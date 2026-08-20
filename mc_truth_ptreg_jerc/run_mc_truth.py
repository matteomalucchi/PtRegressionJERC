#!/usr/bin/env python3
"""Launch the column-based MC truth analysis.

This replaces ``exec.py`` for the MC truth derivation. Differences:

* **no tmux**: the job runs in the foreground and its output goes to the
  terminal (and, with ``--log``, to a file). Use ``nohup``/``screen``/``sbatch``
  yourself if you want it detached;
* **no environment variables to export**: every option is a command-line flag.
  They are written into a run card (a YAML file) inside the output directory,
  and the configuration reads that file;
* **no eta splitting**: one job covers the full eta range, both taggers and
  both the standard and the neutrino-inclusive gen jets.

Examples
--------
Run the full 2023 pre-BPix sample on the PSI T3::

    python run_mc_truth.py -y 2023_preBPix -o /work/$USER/out_mc_truth_2023_preBPix

Quick local test on a few chunks::

    python run_mc_truth.py -y 2023_preBPix -o ./test_out --executor iterative --test

Closure test with the freshly derived corrections::

    python run_mc_truth.py -y 2023_preBPix -o /work/$USER/out_mc_truth_2023_closure --closure
"""

import argparse
import os
import shutil
import subprocess
import sys

import yaml

LOCALDIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_RUNCARD = os.path.join(LOCALDIR, "params", "mc_truth_runcard.yaml")
CONFIG_FILE = os.path.join(LOCALDIR, "mc_truth_columns_config.py")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Run the column-based MC truth analysis",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-y", "--year", default=None, help="Campaign to process (e.g. 2023_preBPix)"
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output directory for the coffea output and the generated run card",
    )
    parser.add_argument(
        "-c",
        "--runcard",
        default=DEFAULT_RUNCARD,
        help="Base run card; command line options are applied on top of it",
    )
    parser.add_argument(
        "--columns-dir",
        default=None,
        help=(
            "Directory (or xrootd URL) for the per-chunk parquet files. "
            "Defaults to the value in the run card; '{year}' is substituted"
        ),
    )
    parser.add_argument(
        "-e",
        "--executor",
        default="dask@T3_CH_PSI",
        help="PocketCoffea executor (e.g. iterative, dask@T3_CH_PSI, condor@lxplus)",
    )
    parser.add_argument(
        "-r",
        "--run-options",
        default=os.path.join(LOCALDIR, "params", "t3_run_options_mc_truth.yaml"),
        help="YAML file with the executor options (cores, memory, chunksize, ...)",
    )
    parser.add_argument(
        "--closure",
        action="store_true",
        help="Apply the derived MC truth corrections on top of the regressions",
    )
    parser.add_argument(
        "--no-pnet", action="store_true", help="Skip the ParticleNet regression"
    )
    parser.add_argument(
        "--no-upart", action="store_true", help="Skip the UParT regression"
    )
    parser.add_argument(
        "--no-neutrinos",
        action="store_true",
        help="Skip the neutrino-inclusive gen jet collection",
    )
    parser.add_argument(
        "--no-reapply-jec",
        action="store_true",
        help="Use the PocketCoffea calibrated jet pT instead of re-applying L2Relative",
    )
    parser.add_argument(
        "--test",
        action="store_true",
        help="Run pocket-coffea in test mode (a couple of chunks, local)",
    )
    parser.add_argument(
        "--process-separately",
        action="store_true",
        help="Process each dataset independently (allows partial recovery)",
    )
    parser.add_argument(
        "--log",
        action="store_true",
        help="Also write the output to <output>/run_mc_truth.log",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Run even if <output>/output_all.coffea already exists",
    )
    parser.add_argument(
        "-n",
        "--dry-run",
        action="store_true",
        help="Only write the run card and print the command",
    )
    return parser.parse_args(argv)


def build_runcard(args):
    """Load the base run card and apply the command-line overrides."""
    with open(args.runcard) as f:
        runcard = yaml.safe_load(f)

    if args.year is not None:
        runcard["year"] = args.year
    if args.columns_dir is not None:
        runcard["columns_output_dir"] = args.columns_dir
    if args.closure:
        runcard["closure"] = True
    if args.no_pnet:
        runcard["pnet"] = False
    if args.no_upart:
        runcard["upart"] = False
    if args.no_neutrinos:
        runcard["neutrinos"] = False
    if args.no_reapply_jec:
        runcard["reapply_jec"] = False

    year = str(runcard["year"])
    if year not in runcard["samples"]:
        raise SystemExit(
            f"Year '{year}' is not defined in the run card. "
            f"Known years: {sorted(runcard['samples'])}"
        )
    # resolve {year} once, so that the workers read an explicit path
    runcard["columns_output_dir"] = runcard["columns_output_dir"].format(year=year)
    return runcard


def write_run_options(args, runcard_path):
    """Copy the executor options next to the output and export the run card.

    Batch executors (condor@lxplus) re-import the configuration file on the
    worker node, so the worker needs to know where the run card is. A single
    ``custom-setup-commands`` entry is enough; the run card must live on a
    shared file system in that case.
    """
    if not args.run_options or not os.path.isfile(args.run_options):
        return None

    dest = os.path.join(args.output, os.path.basename(args.run_options))
    shutil.copy2(args.run_options, dest)

    with open(dest) as f:
        options = yaml.safe_load(f) or {}
    setup_commands = list(options.get("custom-setup-commands", []) or [])
    setup_commands.append(f"export MC_TRUTH_RUNCARD={runcard_path}")
    options["custom-setup-commands"] = setup_commands
    with open(dest, "w") as f:
        yaml.safe_dump(options, f, sort_keys=False)

    return dest


def main(argv=None):
    args = parse_args(argv)
    os.makedirs(args.output, exist_ok=True)

    output_all = os.path.join(args.output, "output_all.coffea")
    if os.path.isfile(output_all) and not args.overwrite:
        print(f"{output_all} already exists, nothing to do (use --overwrite to rerun).")
        return 0

    runcard = build_runcard(args)
    runcard_path = os.path.abspath(os.path.join(args.output, "mc_truth_runcard.yaml"))
    with open(runcard_path, "w") as f:
        yaml.safe_dump(runcard, f, sort_keys=False)
    print(f"Run card written to {runcard_path}")

    run_options = write_run_options(args, runcard_path)

    command = [
        "pocket-coffea",
        "run",
        "--cfg",
        CONFIG_FILE,
        "-o",
        args.output,
    ]
    if args.test:
        command.append("--test")
    else:
        command += ["-e", args.executor]
        if run_options:
            command += ["--custom-run-options", run_options]
    if args.process_separately:
        command.append("--process-separately")

    env = dict(os.environ, MC_TRUTH_RUNCARD=runcard_path)

    print("Running:", " ".join(command))
    if args.dry_run:
        return 0

    if args.log:
        log_path = os.path.join(args.output, "run_mc_truth.log")
        print(f"Logging to {log_path}")
        with open(log_path, "w") as log_file:
            process = subprocess.Popen(
                command, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
            )
            for line in process.stdout:
                line = line.decode(errors="replace")
                sys.stdout.write(line)
                log_file.write(line)
            returncode = process.wait()
    else:
        returncode = subprocess.run(command, env=env).returncode

    if returncode != 0:
        print(f"pocket-coffea failed with exit code {returncode}")
        return returncode

    print("\nDone. Derive the MC truth corrections with:\n")
    print(
        "    python response_plot/response_mc_truth.py "
        f"-i {args.output} -o {os.path.join(args.output, 'mc_truth_plots')}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
