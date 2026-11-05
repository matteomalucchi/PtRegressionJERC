# MC Truth corrections for pT regression

Repository to compute MC Truth corrections and JER MC for regressed pT jets.

This repository is structured as an analysis configurations for [PocketCoffea](https://github.com/PocketCoffea/PocketCoffea/tree/main) and was originally hosted [at this link](https://github.com/matteomalucchi/AnalysisConfigs/tree/legacy-jme-config/configs/jme).

## Overview

This directory contains two independent workflows. They use different PocketCoffea configurations, different launch commands and different post-processing scripts, so they are described in separate sections:

| Workflow | PocketCoffea config | How to launch | Post-processing |
|---|---|---|---|
| [MC Truth corrections](#mc-truth-corrections) (L2Relative) | `cartesian_config.py` | `exec.py` (wraps `pocket-coffea run`) | `response_plot/response.py`, `response_plot/plot_summary_reponse.py` |
| [JER from MC truth](#jer-from-mc-truth) | `jer_mc_config.py` | `pocket-coffea run` directly | `response_plot/plot_jer_mc.py`, `response_plot/plot_jerc_txt.py` |

All the commands below assume that the environment is active (see the main [README](../README.md#activate-the-environment)) and that you are inside the `mc_truth_ptreg_jerc/` directory.

The post-processing scripts are heavy (many cores, a lot of memory, long runtime). Instructions to run them on a batch system (SLURM on the PSI Tier-3, HTCondor on `lxplus`) are given in [Running the post-processing scripts on a batch system](#running-the-post-processing-scripts-on-a-batch-system).

---

## MC Truth corrections

### Running the analysis

To run this over the full dataset for a particular year in each $\eta$ and $p_T$ bin, you can use the following command:

```bash
python exec.py --full [-pnet] [-upart] --dir <dir_name> -y <year> [--lxplus]
```

Where `<dir_name>` is the name of the directory where you want to save the results, and `<year>` is the year you want to run the analysis for. The `--lxplus` flag is used to indicate that you are running this on `lxplus`: jobs are submitted with `-e condor@lxplus`, and `exec.py` creates a temporary run-options file (`params/lxplus_run_options_big_<year>*.tmp.yaml`) that exports the needed environment variables on the condor workers. Without `--lxplus`, jobs are submitted on the PSI Tier-3 with `-e dask@T3_CH_PSI` and `params/t3_run_options_big.yaml`.

Year can be set to:

- 2016_PreVFP
- 2016_PostVFP
- 2017
- 2018
- 2022_preEE
- 2022_postEE
- 2023_preBPix
- 2023_postBPix
- 2024
- 2025

This will save the results in the `dir_name` directory inside the
`output_all.coffea` file.

If running on `lxplus`, there will be an output file for each worker in the `dir_name` directory and you can merge them using:

```bash
cd <dir_name>
pocket-coffea merge-outputs -o output_all.coffea output_job_*.coffea
```

The output file contains 2D histograms for each $\eta$ bin in which the x-axis is the jet $p_T$ response and the y-axis is the jet $p_T$.

### Computing the MC Truth corrections

After running the full dataset, in order to compute the MC Truth corrections, you can use the following command:

```bash
cd response_plot/
python response.py --full -d <dir_name> --histo 
```

To run on SLURM on tier3:

```bash
cd response_plot/
sbatch -p short --account=t3 --time=00:10:00 --mem 15gb --cpus-per-task=32 --wrap="python response.py --full -d  <dir_name> --histo -n 32"
```

To run on HTCondor on `lxplus`, see [HTCondor on lxplus](#htcondor-on-lxplus) (an example for `response.py` is given there).

This will:

- Compute the median of the response in each bin in $\eta$ as a function of $p_T$.
- Get the inverse of the median.
- Fit the inverse of the median with a 6th order polynomial.
- Save the results in the configuration file.

It will also:

- Plot the histograms of the response in each bin in $\eta$ and $p_T$ bin.
- Plot the median of the response in each bin in $\eta$ as a function of $p_T$.
- Plot the inverse of the median in each bin in $\eta$ as a function of $p_T$.
- Plot the resolution of the response in each bin in $\eta$ as a function of $p_T$ using 3 different definitions.

### Closure test

To run the closure test of the corrections you can re-run the analysis with some additional flags:

```bash
python exec.py --full -pnet --dir <dir_name> -y <year> --closure --abs-eta-inclusive [--lxplus]
```

This will run the analysis applying the newly derived corrections which have to be specified in the config file.
Once this is done, you can run the other steps of the analysis to obtain the final plots.

To plot all eta bins on the same plot you can use the following command:

```bash
cd response_plot/
python plot_summary_reponse.py -d <dir_name>
```

This is useful to plot the closure test of the MC Truth corrections in a inclusive way.

---

## JER from MC truth

The JER (Jet Energy Resolution) workflow has three steps:

1. **Run `jer_mc_config.py` with PocketCoffea.** The MC truth corrections (standard, PNet regression, UParT regression, with and without neutrinos) are applied to the jets and, instead of filling histograms, the per-jet columns (response, $p_T$, $\eta$, nPU, $\rho$, ...) are dumped chunk by chunk as parquet files. The `.coffea` files written in the output directory (`-o`) reference these chunks.
2. **Run `response_plot/plot_jer_mc.py`** on the output directory to build the response histograms, fit them, compute the resolution and write the JER txt files (see [JER plotting with `plot_jer_mc.py`](#jer-plotting-with-plot_jer_mcpy)).
3. **Optionally run `response_plot/plot_jerc_txt.py`** on the produced txt files to compare them (see [Plotting JERC txt files with `plot_jerc_txt.py`](#plotting-jerc-txt-files-with-plot_jerc_txtpy)).

`exec.py` is **not** used for this workflow: `pocket-coffea run` is called directly.

### Running `jer_mc_config.py`

#### Before running

- **Year.** The year is read from the `YEAR` environment variable inside `jer_mc_config.py` (default: `2022_preEE`). The allowed values are the same as for the [MC Truth corrections](#running-the-analysis). The year selects the QCD sample and the MC truth correction files that are applied (closure), so the correction files for that year must exist in `params/`.
- **Location of the parquet chunks.** The parquet chunks are **not** written in the `-o` directory, but in the path defined by the `output_chunks_name` variable at the top of `jer_mc_config.py`:

  ```python
  output_chunks_name = f"root://t3dcachedb03.psi.ch:1094//pnfs/psi.ch/cms/trivcat/store/user/mmalucch/out_jer_mc_ptreg/out_jer_closure_pnet_upart_{year}{add_str}/parquet_files"
  ```

  By default this points to the PSI Tier-3 dCache area of the original author. Change it to a location where you have write access before running (see the [lxplus section](#on-lxplus-htcondor) for an EOS example). `add_str` (default `_RecoEta`) is a suffix to distinguish different productions.

#### On the PSI Tier-3 (Dask on SLURM)

```bash
export YEAR=2024
pocket-coffea run \
    --cfg jer_mc_config.py \
    --custom-run-options params/t3_run_options_jer_15Gb.yaml \
    -o /work/<user>/out_jer_mc_ptreg/out_jer_closure_pnet_upart_$YEAR \
    --process-separately \
    -e dask@T3_CH_PSI
```

Key flags:

| Flag | Description |
|---|---|
| `--cfg` | PocketCoffea configuration file (`jer_mc_config.py`) |
| `--custom-run-options` | YAML file with executor settings (cores, memory, walltime, chunk size, etc.) |
| `-o` | Output directory where the per-dataset `.coffea` files are written (the parquet chunks go to `output_chunks_name`) |
| `--process-separately` | Process each dataset independently (enables partial recovery) |
| `-e dask@T3_CH_PSI` | Run on the T3 PSI cluster via Dask (`-e condor@lxplus` on lxplus, `-e iterative` for local testing) |

Since the Dask scheduler runs in your shell for the whole duration of the processing, launch the command inside a `tmux` (or `screen`) session so that it survives a disconnection.

The available run-options files for the JER on the Tier-3 are:

| File | Workers | Memory per worker | Walltime |
|---|---|---|---|
| `params/t3_run_options_jer.yaml` | 30 | 8 GB | 3 h |
| `params/t3_run_options_jer_15Gb.yaml` | 25 | 15 GB | 4 h |
| `params/t3_run_options_jer_15Gb_8h.yaml` | 25 | 15 GB | 8 h |

For example, `params/t3_run_options_jer_15Gb.yaml` contains:

```yaml
scaleout: 25          # number of parallel workers
cores-per-worker: 1
mem-per-worker: "15GB"
disk-per-worker: "1GB"
queue: "standard"
walltime: "04:00:00"
chunksize: 100000
```

Use the 15 GB versions if the workers run out of memory, and the 8 h version if they hit the walltime.

To process several years, run the command once per year, changing `YEAR` (and the `-o` directory) each time:

```bash
for YEAR in 2022_preEE 2022_postEE 2023_preBPix 2023_postBPix 2024; do
    export YEAR
    pocket-coffea run --cfg jer_mc_config.py \
        --custom-run-options params/t3_run_options_jer_15Gb.yaml \
        -o /work/<user>/out_jer_mc_ptreg/out_jer_closure_pnet_upart_$YEAR \
        --process-separately -e dask@T3_CH_PSI
done
```

#### On lxplus (HTCondor)

On `lxplus` the processing is submitted to HTCondor with the `condor@lxplus` executor. Run everything from inside the PocketCoffea apptainer image with the virtual environment active (see the main [README](../README.md#lxplus-installation)).

1. **Change the output location of the parquet chunks.** In `jer_mc_config.py`, set `output_chunks_name` to a directory on EOS where you have write access, e.g.:

   ```python
   output_chunks_name = f"root://eosuser.cern.ch//eos/user/<u>/<user>/out_jer_mc_ptreg/out_jer_closure_pnet_upart_{year}{add_str}/parquet_files"
   ```

   The dCache path at PSI is not writable from the CERN batch nodes.

2. **Create a run-options file that exports `YEAR` on the condor workers.** The condor jobs do not inherit the environment of your shell, so the environment variables read by the configuration must be exported through `custom-setup-commands` (this is what `exec.py --lxplus` does automatically for the MC Truth workflow). Start from `params/lxplus_run_options_big.yaml` (10 workers, 15 GB each, `workday` queue, `local-virtualenv: true`):

   ```bash
   export YEAR=2024
   cp params/lxplus_run_options_big.yaml params/lxplus_run_options_jer_$YEAR.tmp.yaml
   cat >> params/lxplus_run_options_jer_$YEAR.tmp.yaml <<EOT

   custom-setup-commands:
     - export YEAR=$YEAR
   EOT
   ```

   Adjust `scaleout`, `mem-per-worker` and `queue` (the HTCondor `JobFlavour`, see [below](#htcondor-on-lxplus)) to your needs. If the workers do not pick up your editable PocketCoffea installation, add also the `source .../pocket_coffea_env/bin/activate` and `export PYTHONPATH=...` lines that are commented out in `params/lxplus_run_options_big.yaml` to `custom-setup-commands`.

3. **Run** (preferably inside `tmux`, because the submission and the monitoring of the jobs happen in your shell):

   ```bash
   pocket-coffea run \
       --cfg jer_mc_config.py \
       --custom-run-options params/lxplus_run_options_jer_$YEAR.tmp.yaml \
       -o /eos/user/<u>/<user>/out_jer_mc_ptreg/out_jer_closure_pnet_upart_$YEAR \
       --process-separately \
       -e condor@lxplus
   ```

   You can follow the jobs with `condor_q`. If they fail with import errors, check the [Troubleshooting (lxplus)](../README.md#troubleshooting-lxplus) section of the main README.

4. **Outputs.** As for the MC Truth workflow, on lxplus there can be one `.coffea` file per job in the `-o` directory. `plot_jer_mc.py` reads **every** `.coffea` file in its input directory, so you can give it the directory as it is. If you prefer to merge the outputs with `pocket-coffea merge-outputs`, write the merged file to a different directory, otherwise the same events are counted twice.

#### Quick local test

To check that the configuration runs before submitting the full processing, run on a few chunks locally:

```bash
export YEAR=2024
pocket-coffea run --cfg jer_mc_config.py --test -o test_jer_$YEAR
```

Note that the parquet chunks are still written to `output_chunks_name`, so point it to a test location (or a local directory) first.

### JER plotting with `plot_jer_mc.py`

`response_plot/plot_jer_mc.py` is the main plotting script for JER (Jet Energy Resolution) studies from MC truth. It reads the parquet/coffea output produced by the pocket-coffea analysis, builds multi-dimensional histograms of the jet response and auxiliary variables, fits Gaussian distributions to extract the resolution, and produces a full set of publication-quality plots.

#### What the script does

1. **Loads data** — reads all `.coffea` files in the input directory. Each file contains per-chunk arrays of jet variables (response, $p_T$, $\eta$, pile-up, $\rho$, etc.) for each event category and correction type (JEC, Raw, PNet regression, UParT regression, with and without neutrinos).
2. **Builds mapping-variable histograms** — fills N-dimensional histograms for auxiliary variables (`Rho_fixedGridRhoFastjetAll`, reconstructed $p_T$ for each correction) binned in nPU / $\eta$ / gen-$p_T$. A linear fit of $\rho$ vs $\mu$ is performed to convert nPU bin edges into $\rho$ values for the JER text output.
3. **Builds response histograms** — fills per-bin response histograms ($p_T^\text{reco}/p_T^\text{ptcl}$) in a 3-D grid of (nPU, $\eta$, gen-$p_T$). Separate histograms are kept for regular jets and neutrino-inclusive jets.
4. **Fits the response** — fits a Gaussian (with optional tail clipping) to each response histogram to extract the mean and $\sigma$.
5. **Computes JER** — derives the relative resolution $\sigma/\mu$ in each (nPU, $\eta$) slice as a function of gen-$p_T$, and optionally fits the NSC formula: $\sigma/p_T = \sqrt{N^2/p_T^2 + S^2 p_T^d + C^2}$.
6. **Saves results** — writes a coffea file containing all histograms and fit results so they can be reloaded later without re-running the full chain (see `--load` / `--refit`).
7. **Writes JER text files** — exports resolution results in the standard JME text format for each jet correction type and $\eta$ bin.

#### Command-line usage

```bash
cd response_plot/
python plot_jer_mc.py \
    -i <input_dir> \
    -o <output_dir> \
    [-c plot_jer_configs/plot_config_jer_mc_coarsePU_rebinResponse.yaml] \
    [-w <N_workers>] \
    [--histo] \
    [--novars] \
    [-m <max_parquet>] \
    [-l <precomputed.coffea>] \
    [--refit] \
    [--test]
```

#### Arguments

| Argument | Default | Description |
|---|---|---|
| `-i`/`--input-dir` | required | Directory containing `.coffea` output files from pocket-coffea |
| `-o`/`--output` | `./plot_jer_mc/` | Output directory for all plots, coffea file, and log |
| `-c`/`--config` | `plot_jer_configs/plot_config_jer_mc_default.yaml` | Path to the YAML configuration file (see below). When omitted, the built-in default config is used automatically. Any option present in the provided file overwrites the corresponding default; options absent keep their default values. |
| `-w`/`--workers` | `1` | Number of parallel workers for multiprocessing (speeds up histogram filling and Gaussian fitting) |
| `--histo` | off | Also produce 1D/2D histogram plots of the response and mapping variables |
| `--novars` | off | Use the old save format where variation arrays are not stored inside the coffea files |
| `-m`/`--max-parquet` | `None` (all) | Limit the number of parquet files loaded per dataset (useful for quick tests) |
| `-l`/`--load` | `None` | Path to a precomputed coffea file produced by a previous run; skips data loading and histogram filling |
| `--refit` | off | When used with `--load`, re-merge histogram bins according to the current YAML config (coarser binning) and recompute all Gaussian fits before plotting |
| `--test` | off | Run in test mode with reduced bin arrays (overrides the `test` flag in the YAML config) |

#### Typical workflows

**First run from scratch (using the default config):**

```bash
python plot_jer_mc.py -i /work/<user>/out_jer_mc_ptreg/out_jer_closure_pnet_upart_2022_postEE \
    -o plots_2022postEE \
    -w 16 --histo
```

**First run from scratch (with a custom config):**

```bash
python plot_jer_mc.py -i /work/<user>/out_jer_mc_ptreg/out_jer_closure_pnet_upart_2022_postEE \
    -o plots_2022postEE \
    -c plot_jer_configs/plot_config_jer_mc_coarsePU_rebinResponse.yaml \
    -w 16 --histo
```

**Reload and replot with coarser binning (no re-reading of parquet):**

```bash
python plot_jer_mc.py \
    -l plots_2022postEE/plotted_data.coffea \
    -o plots_2022postEE_rebin \
    -c plot_jer_configs/plot_config_jer_mc_coarsePU_rebinResponse.yaml \
    --refit -w 16
```

**Quick test on a subset of files:**

```bash
python plot_jer_mc.py -i /work/<user>/out_jer_... -o test_out \
    -c plot_jer_configs/plot_config_jer_mc_coarsePU_rebinResponse.yaml \
    -m 5
```

(Use `test: true` in the YAML config to also reduce the bin granularity.)

---

### Running the post-processing scripts on a batch system

`plot_jer_mc.py` (and `response.py` for the MC Truth corrections) can take hours and a lot of memory on a full dataset, so they should be run on a batch node rather than interactively. In both cases match the number of workers of the script (`-w` for `plot_jer_mc.py`, `-n` for `response.py`) to the number of requested cores.

#### Example for SLURM on the PSI Tier-3

```bash
cd response_plot/
sbatch -p standard --account=t3 --time=12:00:00 --mem 60gb --cpus-per-task=32 \
    --job-name=jer_mc_$YEAR --output=slurm_jer_mc_%j.out \
    --wrap="python plot_jer_mc.py -i /work/<user>/out_jer_mc_ptreg/out_jer_closure_pnet_upart_$YEAR -o plots_jer_$YEAR -c plot_jer_configs/<config>.yaml -w 32 --histo"
```

#### Example for HTCondor on lxplus

`lxplus` has no `sbatch`; the equivalent is to write a small executable script plus a submit file and submit it with `condor_submit`. The job must run in the same apptainer image and virtual environment used interactively, otherwise the packages are not found. The setup below is generic: the same wrapper is used for `plot_jer_mc.py`, `plot_jerc_txt.py` and `response.py`, only the arguments change.

Put the following files in a directory on **AFS** (HTCondor at CERN cannot write the `log`/`output`/`error` files to EOS), for example `~/condor_jer/`. Inputs and outputs of the scripts themselves can be on EOS.

1. **Wrapper script** `~/condor_jer/run_in_env.sh` (replace the paths with the ones of your installation):

   ```bash
   #!/bin/bash
   # Usage: run_in_env.sh <script.py> [script arguments...]
   # Runs a script of mc_truth_ptreg_jerc/response_plot inside the PocketCoffea environment.
   set -e

   POCKET_COFFEA_DIR=/afs/cern.ch/user/<u>/<user>/PocketCoffea
   RESPONSE_PLOT_DIR=/afs/cern.ch/user/<u>/<user>/PtRegressionJERC/mc_truth_ptreg_jerc/response_plot

   source ${POCKET_COFFEA_DIR}/pocket_coffea_env/bin/activate
   export PYTHONPATH=${POCKET_COFFEA_DIR}:$PYTHONPATH

   cd ${RESPONSE_PLOT_DIR}
   echo "Running on $(hostname) with $(nproc) cores: python $*"
   python "$@"
   ```

   and make it executable:

   ```bash
   chmod +x ~/condor_jer/run_in_env.sh
   ```

2. **Submit file** `~/condor_jer/plot_jer_mc.sub`:

   ```
   universe          = vanilla
   executable        = run_in_env.sh
   arguments         = "plot_jer_mc.py -i /eos/user/<u>/<user>/out_jer_mc_ptreg/out_jer_closure_pnet_upart_$(YEAR) -o /eos/user/<u>/<user>/plots_jer_mc/plots_jer_$(YEAR) -c plot_jer_configs/<config>.yaml -w 16 --histo"

   # Run inside the PocketCoffea image (same as the interactive apptainer shell)
   +SingularityImage = "/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-analysis/general/pocketcoffea:lxplus-el9-stable"

   request_cpus      = 16
   request_memory    = 32000MB
   +JobFlavour       = "workday"

   log               = logs/plot_jer_mc_$(YEAR).$(ClusterId).log
   output            = logs/plot_jer_mc_$(YEAR).$(ClusterId).out
   error             = logs/plot_jer_mc_$(YEAR).$(ClusterId).err

   YEAR = 2024
   queue
   ```

3. **Submit** from `lxplus` (outside or inside the apptainer shell):

   ```bash
   cd ~/condor_jer
   mkdir -p logs
   condor_submit plot_jer_mc.sub
   ```

   To override the year without editing the file:

   ```bash
   condor_submit YEAR=2022_postEE plot_jer_mc.sub
   ```

   To submit one job per year in a single command, replace the last two lines of the submit file (`YEAR = 2024` and `queue`) with:

   ```
   queue YEAR in (2022_preEE 2022_postEE 2023_preBPix 2023_postBPix 2024)
   ```

4. **Monitor and manage** the jobs:

   ```bash
   condor_q
   ```

   ```bash
   condor_tail -f <job_id>
   ```

   ```bash
   condor_rm <job_id>
   ```

   ```bash
   tail -f logs/plot_jer_mc_2024.*.out
   ```

Notes:

- **Runtime.** `+JobFlavour` sets the maximum wall time: `espresso` (20 min), `microcentury` (1 h), `longlunch` (2 h), `workday` (8 h), `tomorrow` (1 day), `testmatch` (3 days), `nextweek` (1 week). Use `tomorrow` or longer for a full first run from scratch; `--load`/`--refit` re-runs are much faster.
- **Memory.** On the CERN batch system about 2 GB of memory are guaranteed per requested core; if you need more memory, increase `request_cpus` (and `-w`) together with `request_memory`.
- **Reload and refit on condor.** Just change `arguments`, e.g.

  ```
  arguments = "plot_jer_mc.py -l /eos/user/<u>/<user>/plots_jer_mc/plots_jer_$(YEAR)/plotted_data.coffea -o /eos/user/<u>/<user>/plots_jer_mc/plots_jer_$(YEAR)_rebin -c plot_jer_configs/<config>.yaml --refit -w 16"
  ```

- **`response.py` (MC Truth corrections).** Use the same wrapper with a different submit file, e.g. `~/condor_jer/response.sub`:

  ```
  universe          = vanilla
  executable        = run_in_env.sh
  arguments         = "response.py --full -d <dir_name> --histo -n 32"
  +SingularityImage = "/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-analysis/general/pocketcoffea:lxplus-el9-stable"
  request_cpus      = 32
  request_memory    = 64000MB
  +JobFlavour       = "longlunch"
  log               = logs/response.$(ClusterId).log
  output            = logs/response.$(ClusterId).out
  error             = logs/response.$(ClusterId).err
  queue
  ```

  `<dir_name>` must be an absolute path. The minimal `response_plot/job.submit` in the repository calls `/usr/bin/python3` directly, without the apptainer image and the virtual environment, so the PocketCoffea packages are not available there; prefer the wrapper above.
- **Import errors on the workers** (e.g. `ModuleNotFoundError: No module named 'utils_configs'`): see [Troubleshooting (lxplus)](../README.md#troubleshooting-lxplus) in the main README.

---

### Plotting JERC txt files with `plot_jerc_txt.py`

`response_plot/plot_jerc_txt.py` reads one or more JME-format resolution (or L2Relative correction) txt files and produces publication-quality plots of the resolution (or jet response) as a function of jet $\eta$, evaluated at a set of user-chosen $p_T$ values.

#### What the script does

1. **Parses the JERC txt format** — reads the file header (bin variables, x variable, formula string) and each row of bin edges and formula parameters. Rows with non-finite parameters are logged and skipped.
2. **Evaluates the formula** — for every combination of (JetPt, $\eta$ bin, $\rho$ bin), the ROOT `TFormula` string is evaluated numerically via numpy. Unphysical combinations (where $E_\text{jet} = p_T \cosh|\eta| \geq E_\text{beam}/2$) are silently dropped.
3. **Produces resolution-vs-$\eta$ plots** — one plot per input file, with:
   - **Color + marker** per $p_T$ evaluation value.
   - **Alpha shading** per $\rho$ bin (lighter = lower $\rho$, darker = higher $\rho$). When no $\rho$ binning is present a single opaque marker per $p_T$ is used.
   - **CMS detector-region lines and labels** (Barrel, EC1, EC2, HF) overlaid on the $\eta$ axis.
   - A full-range plot (all $\eta$) and a restricted-range plot ($|\eta| < $ `--eta-max`, default 2.5) are saved for each file.
4. **Handles L2Relative correction files** — when the txt header declares `Correction L2Relative`, the script inverts the correction ($1/c$) and plots it as the simulated jet response with a horizontal reference line at 1.
5. **Produces pairwise ratio plots** — for every pair of input files, plots the point-wise ratio (numerator / denominator) as a function of $\eta$ with a horizontal line at 1. Ratio plots are also saved in the restricted-$\eta$ range. Use `--no-ratio` to skip this step.

#### Command-line usage

```bash
cd response_plot/

# single file
python plot_jerc_txt.py Run3Summer22_V1_NSC_MC_PtResolution_AK4PFPuppi.txt \
    --pt-values 80 150 300 600 -o plots/

# all txt files in a directory (shell expands the glob)
python plot_jerc_txt.py /path/to/txt/*.txt --pt-values 80 150 300 600 -o plots/
```

#### Arguments

| Argument | Default | Description |
|---|---|---|
| `input` | required | One or more JERC resolution/correction txt files; shell globs (`*.txt`) are accepted |
| `--pt-values` | `50 100 300 600 1000 3000` | JetPt [GeV] values at which to evaluate the resolution formula |
| `-o`/`--output` | `.` | Output directory for all plots |
| `--formats` | `png pdf svg` | Output file formats |
| `--no-ratio` | off | Skip pairwise ratio plots |
| `--eta-max` | `2.5` | $|\eta|$ upper limit for the restricted-range plots |

#### Output files

For each input file `<stem>.txt` the script writes:

- `<output>/<stem>.<fmt>` — full-$\eta$ resolution (or response) plot.
- `<output>/<stem>_abseta_lt_<eta_max>.<fmt>` — restricted-$\eta$ version.

For each pair of files `<stem_a>` and `<stem_b>`:

- `<output>/ratio_<stem_a>__over__<stem_b>.<fmt>` — full-$\eta$ ratio plot.
- `<output>/ratio_<stem_a>__over__<stem_b>_abseta_lt_<eta_max>.<fmt>` — restricted-$\eta$ ratio.

---

### Plot configuration file (`plot_jer_configs/`)

The YAML configuration file controls every aspect of the binning, variables, response histograms, and plot style. Several pre-made configurations exist in `response_plot/plot_jer_configs/`.

The file `plot_jer_configs/plot_config_jer_mc_default.yaml` is the built-in default loaded automatically when `-c` is not specified. Every option in that file is annotated with an inline comment explaining its purpose. Pass a custom config with `-c` to override only the options you need; the rest will keep their default values.

#### Boolean flags

| Key | Default | Description |
|---|---|---|
| `test` | `false` | Replace all bin arrays with the reduced `test_*` arrays for fast iteration |
| `mixed_mode` | `true` | Process regular and neutrino-inclusive response variables together in a single `mixed` mode; when `false` they are processed as two separate `regular` and `neutrino` modes |
| `resolution_vs_pt_gen` | `false` | Plot JER as a function of gen-$p_T$ |
| `resolution_vs_pt_reco` | `false` | Plot JER as a function of mean reco-$p_T$ (mapped x-axis) |
| `histograms_map` | `false` | Save/plot the mapping-variable histograms |
| `histograms_response` | `false` | Save/plot the 1D response histograms per bin |
| `gaussian_fit_resolution` | `true` | Derive resolution from a Gaussian fit to the response histogram |
| `gaussian_fit_cut_tails` | `[0.87]` | List of CI fractions to try when trimming the fit window around the peak; set to `false` to disable. The CI that yields the best p-value is kept per bin |
| `gaussian_fit_max_sigma_rel_err` | `1.0` | Reject bins where the fitted `sigma_err / sigma` exceeds this value (resolution becomes NaN) |
| `min_events_for_fit` | `50` | Minimum number of events in a 1D response slice before any fit or CI estimation is attempted |
| `rebin_max_factor` | `None` | Maximum rebin factor for the auto-rebinning algorithm. Prevents low-statistics slices (e.g. high pT, forward η) from being rebinned so coarsely that the Gaussian fit over-estimates the resolution. Set to `null` to disable (original behaviour). Typical range: 8–12 |
| `ci_conf_level` | `0.87` | Confidence level used by the CI estimator when `gaussian_fit_resolution: false` |
| `save_response_summary` | `false` | Save a JSON summary of response histograms and fit results for later use in the JME text-file output |

#### Fitting backend and convergence options

All fits (per-bin Gaussian resolution, NSC resolution-vs-$p_T$, linear $\rho$-vs-$\mu$ mapping) share a single backend selected by `fit_backend`. The default `curve_fit` keeps the previous behaviour; switch to `iminuit` for the MINUIT minimizer ([documentation](https://scikit-hep.org/iminuit/)), which is often more robust for the NSC and Gaussian fits.

```yaml
fit_backend: curve_fit   # "curve_fit" (default) or "iminuit"
```

**`iminuit_options`** (used only when `fit_backend: iminuit`):

| Sub-key | Default | Description |
|---|---|---|
| `optimizer` | `migrad` | Minimiser to run: `"migrad"` (recommended), `"simplex"`, or `"scan"` |
| `run_simplex_first` | `false` | Run a SIMPLEX pass before MIGRAD to escape shallow local minima |
| `strategy` | `1` | MINUIT strategy level: `0` = fast, `1` = default, `2` = accurate (slower, better error estimates) |
| `tol` | `null` | EDM convergence tolerance; `null` uses the iminuit default |
| `max_calls` | `0` | Maximum function calls per minimisation step; `0` = no limit |
| `use_minos` | `false` | Run MINOS after MIGRAD to compute asymmetric uncertainties |

**`curve_fit_options`** (used only when `fit_backend: curve_fit`):

| Sub-key | Default | Description |
|---|---|---|
| `method` | `null` | scipy method: `null` (auto: LM when unbounded, TRF when bounded), `"lm"`, `"trf"`, or `"dogbox"` |
| `maxfev` | `100000` | Maximum number of function evaluations |

#### NSC fit options

Options controlling the NSC formula ($\sigma/p_T = \sqrt{N^2/p_T^2 + S^2 p_T^d + C^2}$) fit to the resolution-vs-$p_T$ curves:

| Key | Default | Description |
|---|---|---|
| `nsc_fit_parameters` | `[3.5, 0.5, 0.05, -1.0]` | Initial values `[N, S, C, d]` for the NSC fit |
| `nsc_fit_parameters_by_eta` | `null` | List of `{eta_range: [lo, hi], nsc_fit_parameters: [...]}` rules; per-eta override for initial parameters |
| `nsc_fit_parameter_limits` | all `null` | Lower and upper bounds per parameter; `null` entries are ±∞. Set the whole key to `null` to disable all limits |
| `nsc_fit_parameter_limits_by_eta` | `null` | List of `{eta_range: [lo, hi], nsc_fit_parameter_limits: {lower: [...], upper: [...]}}` rules; per-eta override for parameter limits |
| `nsc_fix_parameters` | all `null` | Per-parameter list: float = fix to that value, `null` = free. Example: `[null, null, null, -1.0]` fixes `d=-1.0` |
| `nsc_fix_parameters_by_eta` | `null` | List of `{eta_range: [lo, hi], nsc_fix_parameters: [...]}` rules; per-eta override for fixed parameters |
| `nsc_fit_x_clip` | `null` | `[min, max]` range to clip the x-axis before fitting (e.g. `[15, 3000]`); `null` disables clipping |
| `nsc_fit_x_clip_by_eta` | `null` | List of `{eta_range: [lo, hi], x_clip: [min, max]}` rules; the first rule whose `eta_range` contains the current bin's \|η\| midpoint wins and overwrites `nsc_fit_x_clip` for that bin |
| `plot_fit_diagnostics` | `false` | When `true`, draw a black dotted line showing the NSC formula at the initial parameters (p₀), added to the legend, for the **first response variable only** |

These keys can also be set **per response variable** to override the global values for that variable only.

`nsc_fit_x_clip_by_eta` example — different x ranges for barrel, endcap, and forward:

```yaml
# global (applies to all response variables)
nsc_fit_x_clip_by_eta:
  - eta_range: [0.0, 1.3]
    x_clip: [15, 3000]
  - eta_range: [1.3, 2.5]
    x_clip: [20, 3000]
  - eta_range: [2.5, 5.2]
    x_clip: [30, 3000]

# per response variable (takes priority over the global list)
response_variables:
  MatchedJets_ResponsePNet:
    nsc_fit_x_clip_by_eta:
      - eta_range: [2.5, 5.2]
        x_clip: [40, 3000]
```

Rules are matched by the absolute-η midpoint of the bin; the first matching rule wins.  Bins whose midpoint falls in no rule fall back to `nsc_fit_x_clip`.

#### Bin arrays

```yaml
pu_bins: [0, 20, 30, 40, 50, 60, 100]   # nPU bin edges
jet_eta_bins: [-5.191, ..., 5.191]       # η bin edges (JME standard)
jet_pt_bins: [8, 10, ..., 5000]          # gen-pT bin edges [GeV]
```

Reduced arrays used when `test: true`:

```yaml
test_pu_bins: [0, 50, 100]
test_jet_eta_bins: [-5.191, 0.0, 5.191]
test_jet_pt_bins: [8, 50, 120, ...]
```

#### Response histogram settings

```yaml
n_bins_response: 150          # number of bins in the response (pT_reco/pT_gen) histogram
bin_limits_response: [0.5, 1.5]  # [min, max] of the response axis
y_lim_resolution: [0, 0.35]   # y-axis range for resolution plots
```

#### `year_map`

Maps the year tag extracted from dataset names to the campaign string used in JME text-file names:

```yaml
year_map:
  "2022_preEE": Summer22
  "2022_postEE": Summer22EE
  "2023_preBPix": Summer23
  ...
```

#### `bin_variables` and `bin_variables_neutrino`

Define the variables used as binning axes for the response histograms. Each entry is a column name from the coffea arrays:

```yaml
bin_variables:
  Pileup_nPU:
    event_var: true          # true → one value per event; false → one value per jet
    bin_edges: *pu_bins      # YAML anchor referencing the bin array above
    label: '$\mu$'           # axis label for plots
    name_plot: nPU           # short name used to build plot filenames
    resolution_x_variable: false  # true → this is the x-axis of the JER vs pT plot
    txt_map_to: Rho_fixedGridRhoFastjetAll  # variable to map this axis to in txt output
  MatchedJets_eta:
    event_var: false
    bin_edges: *jet_eta_bins
    label: '$\eta$'
    name_plot: jet_eta
    resolution_x_variable: false
    txt_name: JetEta         # name used in the JME text-file header
  MatchedJets_pt:
    event_var: false
    bin_edges: *jet_pt_bins
    label: '$p_{\mathrm{T}}^{ptcl}$ [GeV]'
    name_plot: jet_gen_pt
    resolution_x_variable: true   # this axis becomes the x-axis of the JER plot
    txt_name: JetPt
```

`bin_variables_neutrino` mirrors `bin_variables` but for the neutrino-inclusive jet collection (`MatchedJetsNeutrino_*`).

#### `mapping_variables`

Variables that are profiled and fit (e.g. $\rho$ vs $\mu$, reco-$p_T$ vs gen-$p_T$). Used to map bin centres to mean values for the x-axis of plots and for the JME text output:

```yaml
mapping_variables:
  Rho_fixedGridRhoFastjetAll:
    event_var: true
    label: '$\rho$ [GeV/area]'
    name_plot: rho
    N_bins: 1000
    bin_vars: [Pileup_nPU]    # bins in which the profile is computed
    linear_fit: true          # fit a linear model ρ = m·μ + b per η/pT bin
    txt_name: Rho
    rebin_for_plotting: true
  MatchedJets_JetPtJEC:
    label: '$p_{\mathrm{T}}^{\mathrm{reco}}$ [GeV]'
    name_plot: jet_reco_pt
    N_bins: 1000
    bin_limits: [0, 6000]
    bin_vars: [Pileup_nPU, MatchedJets_eta, MatchedJets_pt]
    legend_name: JEC
    rebin_for_plotting: true
```

Add or remove entries to include additional correction types or auxiliary variables.

#### `response_variables`, `response_variables_neutrino`, `response_variables_mixed`

Define which response columns to read and how to plot them. Each entry corresponds to one correction type.

The basic plotting / labelling keys are:

| Key | Description |
|---|---|
| `label` | Axis label used on response histograms |
| `name_plot` | Short name used in plot file names and for the `plot_settings` lookup |
| `N_bins`, `bin_limits` | Number of bins and `[lo, hi]` range of the response axis |
| `bin_vars` | Ordered list of binning variables (must be keys of `bin_variables` / `bin_variables_neutrino`) |
| `legend_name` | Legend entry in overlay plots |
| `map_x_variable` | Mapping variable whose per-bin mean is used as the x-axis of the resolution-vs-reco-$p_T$ plot |
| `is_reference` | If `true` this curve becomes the denominator of the ratio panel |
| `txt_jet_name` | Jet-algorithm tag inserted in the JME text-file name (e.g. `AK4PFPuppi`) |
| `rebin_for_plotting` | If `true` the stored fine-binned response is rebinned to roughly `min_events_per_bin` entries per bin before fitting/plotting |
| `eta_max` | Optional upper limit on $|\eta|$ for this entry; bins beyond it are skipped (used for AK8 / UParT in the high-$\eta$ region) |

In addition, each response variable supports a set of uncertainty- and
mean-computation controls that override the defaults computed inside
`_compute_resolution_from_histogram`. They are applied per bin in the
order listed below. All default to off so existing configs keep their
previous behaviour.

| Key | Type | Default | Effect |
|---|---|---|---|
| `normalize_by_response_mean` | bool | `false` | Plot the relative resolution `sigma / <R>` instead of `sigma`. The uncertainty is propagated as `(sigma / <R>) * sqrt((sigma_err / sigma)^2 + (mean_err / <R>)^2)` |
| `add_min_err` | float or string | `"0.0"` | Add a minimal uncertainty floor in quadrature: `new_err = sqrt(err^2 + add_min_err^2)`. For example `"0.0005"` enforces a 0.05% floor |
| `add_chi2_ndof` | bool | `false` | Multiply the uncertainty by `sqrt(chi2 / ndf)` wherever the per-bin Gaussian fit's reduced chi-square exceeds 1. Requires `gaussian_fit_resolution: true` |
| `additional_uncertainty` | float or string | `"0.0"` | Add a relative uncertainty in quadrature: `new_err = sqrt(err^2 + (v * sigma)^2)` |
| `response_mean_method` | string | `"auto"` | Selects how the mean of the response is computed per bin. Affects the value used by `normalize_by_response_mean`. Allowed values: `"auto"` (Gaussian fit when the resolution itself comes from a successful fit, binned moments otherwise), `"gaussian_fit"` (force the Gaussian fit mean), `"histogram"` (first moment of the binned 1D response slice), `"mean_storage"` (mean read from a `hist.storage.Mean()` histogram filled with the un-binned response values; built automatically by `create_ND_histo`) |

Example with all options shown:

```yaml
response_variables:
  MatchedJets_ResponseJEC:
    label: '$p_{\mathrm{T}}^{\mathrm{reco}}$/$p_{\mathrm{T}}^{ptcl}$'
    name_plot: response_jec
    N_bins: *n_bins_response
    bin_limits: *bin_limits_response
    bin_vars: [Pileup_nPU, MatchedJets_eta, MatchedJets_pt]
    legend_name: JEC                       # label in the legend
    map_x_variable: MatchedJets_JetPtJEC   # mapping variable for the reco-pT x-axis
    is_reference: true                     # reference curve in overlay/ratio plots
    txt_jet_name: AK4PFPuppi               # jet-algorithm tag in the JME text-file name
    rebin_for_plotting: true               # rebin the fine-binned response before plotting
    eta_max: 5.0                           # optional |eta| cut for this entry
    normalize_by_response_mean: true       # plot sigma / <R> instead of bare sigma
    additional_uncertainty: "0.2"          # add (0.2 * sigma) in quadrature to the error
    add_min_err: "0.0005"                  # minimal uncertainty floor of 0.05%
    add_chi2_ndof: true                    # multiply error by sqrt(chi2/ndf) where >1
    response_mean_method: "mean_storage"   # auto | gaussian_fit | histogram | mean_storage
```

`response_variables_mixed` combines regular and neutrino-inclusive entries and is used when `mixed_mode: true`.

#### `plot_settings`

Controls the color and marker style for each correction type in overlay plots. Keys are matched against the suffix of the response variable name:

```yaml
plot_settings:
  JEC:              { color: darkorange, fmt: "o" }
  Raw:              { color: green,      fmt: "s" }
  PNetReg:          { color: darkred,    fmt: "<" }
  UparTReg:         { color: tomato,     fmt: "^" }
  PNetRegNeutrino:  { color: darkblue,   fmt: ">" }
  UparTRegNeutrino: { color: cornflowerblue, fmt: "v" }
```
