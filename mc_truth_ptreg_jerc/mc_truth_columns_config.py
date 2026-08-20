"""PocketCoffea configuration for the column-based MC truth derivation.

Run it with::

    python run_mc_truth.py -y 2023_preBPix -o /work/<user>/out_mc_truth_2023_preBPix

or directly with::

    pocket-coffea run --cfg mc_truth_columns_config.py -o <outdir> \
        -e dask@T3_CH_PSI --custom-run-options params/t3_run_options_mc_truth.yaml

Everything that used to be steered through environment variables (``YEAR``,
``PNET``, ``UPART``, ``CLOSURE``, ``SIGN``, ``CARTESIAN``, ``FLAVSPLIT``, ...)
now lives in a single YAML run card, ``params/mc_truth_runcard.yaml``. The only
environment variable still understood is ``MC_TRUTH_RUNCARD``, which
``run_mc_truth.py`` sets automatically to point at the run card it generated
from the command line; you never have to export it yourself.

Compared to ``mc_truth_config.py`` / ``cartesian_config.py`` this configuration:

* runs ParticleNet **and** UParT in the same job;
* dumps flat columns into parquet files instead of filling one histogram per
  eta/pT category, so the eta splitting (``neg1``, ``neg2``, ...) is gone;
* keeps the gen jet flavour so that the flavour splitting is done at plotting
  time by ``response_plot/response_mc_truth.py``.
"""

import os
import yaml

from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.parameters.cuts import passthrough
from pocket_coffea.lib.columns_manager import ColOut
from pocket_coffea.lib.cut_definition import Cut
from pocket_coffea.parameters import defaults

import mc_truth_ptreg_jerc.custom_cut_functions as cut_functions
from mc_truth_ptreg_jerc.mc_truth_workflow import MCTruthProcessor
from mc_truth_ptreg_jerc.mc_truth_corrections import load_l2relative_txt
from mc_truth_ptreg_jerc.params.l2relative_txt import get_corr_file

localdir = os.path.dirname(os.path.abspath(__file__))

# -----------------------------------------------------------------------------
# Run card
# -----------------------------------------------------------------------------
runcard_path = os.environ.get(
    "MC_TRUTH_RUNCARD", os.path.join(localdir, "params", "mc_truth_runcard.yaml")
)
with open(runcard_path) as f:
    runcard = yaml.safe_load(f)

print(f"[mc_truth_columns_config] run card: {runcard_path}")

year = str(runcard["year"])
sample = runcard["samples"][year]
columns_output_dir = runcard["columns_output_dir"].format(year=year)

print(f"[mc_truth_columns_config] year:    {year}")
print(f"[mc_truth_columns_config] sample:  {sample}")
print(f"[mc_truth_columns_config] columns: {columns_output_dir}")

# -----------------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------------
default_parameters = defaults.get_default_parameters()
defaults.register_configuration_dir("config_dir", localdir + "/params")

parameters = defaults.merge_parameters_from_files(
    default_parameters,
    f"{localdir}/params/object_preselection.yaml",
    f"{localdir}/params/jets_calibration.yaml",
    update=True,
)

# -----------------------------------------------------------------------------
# MC truth corrections
# -----------------------------------------------------------------------------
(
    corr_files_pnetreg,
    corr_files_upartreg,
    corr_files_pnetreg_neutrino,
    corr_files_upartreg_neutrino,
    corr_files_jec,
) = get_corr_file(os.path.join(localdir, "params"))

mc_truth_corr = None
if runcard.get("reapply_jec", True):
    print(f"[mc_truth_columns_config] re-applying JEC {corr_files_jec[year]}")
    mc_truth_corr = load_l2relative_txt(corr_files_jec[year])

closure_corrections = {}
if runcard.get("closure", False):
    closure_files = {
        "PNetReg": corr_files_pnetreg[year],
        "PNetRegNeutrino": corr_files_pnetreg_neutrino[year],
        "UparTReg": corr_files_upartreg[year],
        "UparTRegNeutrino": corr_files_upartreg_neutrino[year],
    }
    for name, path in closure_files.items():
        print(f"[mc_truth_columns_config] closure test {name}: {path}")
        closure_corrections[name] = load_l2relative_txt(path)

# -----------------------------------------------------------------------------
# Columns dumped to parquet
# -----------------------------------------------------------------------------
# Kinematics + flavour of the matched gen jet, and the reco pT / response of
# every correction. `flatten=False` keeps the jagged (per-event) structure, so
# that event level variables can be broadcast to the jets when reading back.
matched_jets_columns = [
    "pt",  # gen jet pT
    "eta",  # gen jet eta
    "RecoEta",
    "RecoPhi",
    "partonFlavour",
    "hadronFlavour",
    "ResponseJEC",
    "JetPtJEC",
    "ResponseRaw",
    "JetPtRaw",
]
matched_jets_neutrino_columns = list(matched_jets_columns)

if runcard.get("pnet", True):
    matched_jets_columns += ["ResponsePNetReg", "JetPtPNetReg"]
    matched_jets_neutrino_columns += [
        "ResponsePNetRegNeutrino",
        "JetPtPNetRegNeutrino",
    ]
if runcard.get("upart", True):
    matched_jets_columns += ["ResponseUparTReg", "JetPtUparTReg"]
    matched_jets_neutrino_columns += [
        "ResponseUparTRegNeutrino",
        "JetPtUparTRegNeutrino",
    ]

columns = [
    ColOut("Rho", ["fixedGridRhoFastjetAll"]),
    ColOut("Pileup", ["nPU"]),
    ColOut("MatchedJets", matched_jets_columns, flatten=False),
]
if runcard.get("neutrinos", True):
    columns.append(
        ColOut("MatchedJetsNeutrino", matched_jets_neutrino_columns, flatten=False)
    )

# -----------------------------------------------------------------------------
# Preselection: reconstructed and generated primary vertex close to each other
# -----------------------------------------------------------------------------
PV_presel = Cut(
    name="PV_presel",
    params={"distance": runcard.get("pv_presel_distance", 0.2)},
    function=cut_functions.PV_presel_cuts,
)

# -----------------------------------------------------------------------------
# Configurator
# -----------------------------------------------------------------------------
cfg = Configurator(
    parameters=parameters,
    datasets={
        "jsons": [
            os.path.join(localdir, json_file)
            for json_file in runcard["dataset_jsons"]
        ],
        "filter": {
            "samples": [sample],
            "samples_exclude": [],
            "year": [year],
        },
        "subsamples": {},
    },
    workflow=MCTruthProcessor,
    workflow_options={
        "donotscale_sumgenweights": True,
        "dump_columns_as_arrays_per_chunk": columns_output_dir,
        "DeltaR_matching": runcard.get("deltaR_matching", 0.2),
        "GenJetPtCut": runcard.get("genjet_pt_cut", 0),
        "n_jets": runcard.get("n_jets", 3),
        "com_energy": runcard.get("com_energy", 13.6),
        "pnet": runcard.get("pnet", True),
        "upart": runcard.get("upart", True),
        "neutrinos": runcard.get("neutrinos", True),
        "mc_truth_corr": mc_truth_corr,
        "closure_corrections": closure_corrections,
        "regression_value_when_invalid": runcard.get(
            "regression_value_when_invalid", "raw"
        ),
        "invalid_regression_eta": runcard.get("invalid_regression_eta", 4.7),
    },
    calibrators=[],
    skim=[],
    preselections=[PV_presel],
    categories={
        "baseline": [passthrough],
    },
    weights={
        "common": {"inclusive": [], "bycategory": {}},
        "bysample": {},
    },
    variations={
        "weights": {
            "common": {"inclusive": [], "bycategory": {}},
            "bysample": {},
        }
    },
    variables={},
    columns={
        "common": {
            "inclusive": columns,
        },
    },
)
