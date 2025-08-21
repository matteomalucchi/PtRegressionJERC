import os

from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.parameters.cuts import passthrough
from pocket_coffea.lib.columns_manager import ColOut
from pocket_coffea.lib.cut_functions import get_HLTsel
from pocket_coffea.lib.categorization import CartesianSelection, MultiCut
from pocket_coffea.parameters.histograms import (
    met_hists,
    muon_hists,
    count_hist,
    HistConf,
    Axis,
)
from pocket_coffea.parameters import defaults
from pocket_coffea.lib.calibrators.common.common import (
    JetsCalibrator,
    JetsPtRegressionCalibrator,
)

from workflow_Calibrator import METProcessor
from configs.jme.cuts import PV_presel
from custom_cuts import dimuon_presel, at_least_one_jet

localdir = os.path.dirname(os.path.abspath(__file__))

# Loading default parameters

default_parameters = defaults.get_default_parameters()
defaults.register_configuration_dir("config_dir", localdir + "/params")

year = os.environ.get("YEAR", "2022_preEE")
# default_parameters.lepton_scale_factors.electron_sf["apply_ele_scale_and_smearing"][year] = False


# adding object preselection
parameters = defaults.merge_parameters_from_files(
    default_parameters,
    f"{localdir}/params/object_preselection.yaml",
    f"{localdir}/params/triggers.yaml",
    f"{localdir}/params/jets_calibration_Calibrator.yaml",
    update=True,
)


common_cats = {
    "baseline": [passthrough],
}

cfg = Configurator(
    parameters=parameters,
    datasets={
        "jsons": [
            # f"{localdir}/datasets/QCD.json",
            f"{localdir}/datasets/QCD_PNetReg15.json",
            f"{localdir}/datasets/DYJetsToLL_M-50_redirector.json",
        ],
        "filter": {
            "samples": [
                (
                    "DYJetsToLL_M-50"
                    # "DYJetsToLL_M-50_local"
                )
            ],
            "samples_exclude": [],
            # "year": [year],
        },
        "subsamples": {},
    },
    workflow=METProcessor,
    workflow_options={
        # "donotscale_sumgenweights": True,
        "only_physical_jet": True,
        "rescale_MET_with_regressed_pT": True,
        "jec_pt_threshold": 15.0,
    },
    skim=[
        get_HLTsel(primaryDatasets=["SingleMuon"]),
    ],
    preselections=[
        PV_presel,
        dimuon_presel,
        at_least_one_jet,
        # get_nObj_min(2, coll="MuonGood")
    ],
    categories={
        **common_cats,
    },
    weights={
        "common": {
            "inclusive": [
                "genWeight",
                "lumi",
                "XS",
            ],
            "bycategory": {},
        },
        "bysample": {},
    },
    calibrators=[JetsPtRegressionCalibrator, JetsCalibrator],
    variations={
        "weights": {
            "common": {
                "inclusive": [],
                "bycategory": {},
            },
            "bysample": {},
        }
    },
    variables={},
    columns={
        "common": {
            "inclusive": [
                ColOut("ll", ["mass", "pt", "eta", "phi"]),
                ColOut("GenMET", ["pt", "phi"]),
                # ColOut("GenMETPlusNeutrino", ["pt", "phi"]),
                ColOut("RawPuppiMET", ["pt", "phi"]),
                ColOut("RawPuppiMETType1", ["pt", "phi"]),
                ColOut("RawPuppiMETPNet", ["pt", "phi"]),
                ColOut("RawPuppiMETPNetPlusNeutrino", ["pt", "phi"]),
                ColOut("PuppiMET", ["pt", "phi"]),
                ColOut("PuppiMETType1", ["pt", "phi"]),
                ColOut("PuppiMETPNet", ["pt", "phi"]),
                ColOut("PuppiMETPNetPlusNeutrino", ["pt", "phi"]),
                ColOut(
                    "RawPuppiMET_MuonGood",
                    ["pt", "phi", "u_perp_predict", "u_paral_predict", "response"],
                ),
                ColOut(
                    "RawPuppiMETType1_MuonGood",
                    ["pt", "phi", "u_perp_predict", "u_paral_predict", "response"],
                ),
                ColOut(
                    "RawPuppiMETPNet_MuonGood",
                    ["pt", "phi", "u_perp_predict", "u_paral_predict", "response"],
                ),
                ColOut(
                    "RawPuppiMETPNetPlusNeutrino_MuonGood",
                    ["pt", "phi", "u_perp_predict", "u_paral_predict", "response"],
                ),
                ColOut(
                    "PuppiMET_MuonGood",
                    ["pt", "phi", "u_perp_predict", "u_paral_predict", "response"],
                ),
                ColOut(
                    "PuppiMETType1_MuonGood",
                    ["pt", "phi", "u_perp_predict", "u_paral_predict", "response"],
                ),
                ColOut(
                    "PuppiMETPNet_MuonGood",
                    ["pt", "phi", "u_perp_predict", "u_paral_predict", "response"],
                ),
                ColOut(
                    "PuppiMETPNetPlusNeutrino_MuonGood",
                    ["pt", "phi", "u_perp_predict", "u_paral_predict", "response"],
                ),
            ]
        },
        "bysample": {},
    },
)
