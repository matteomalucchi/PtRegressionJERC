import os

from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.parameters.cuts import passthrough
from pocket_coffea.lib.columns_manager import ColOut
from pocket_coffea.lib.categorization import CartesianSelection, MultiCut
from pocket_coffea.parameters.histograms import met_hists
from pocket_coffea.parameters import defaults

from workflow import METProcessor
from configs.jme.cuts import PV_presel


localdir = os.path.dirname(os.path.abspath(__file__))

# Loading default parameters

default_parameters = defaults.get_default_parameters()
defaults.register_configuration_dir("config_dir", localdir + "/params")


year = os.environ.get("YEAR", "2023_preBPix")
# adding object preselection
parameters = defaults.merge_parameters_from_files(
    default_parameters,
    f"{localdir}/params/object_preselection.yaml",
    f"{localdir}/params/jets_calibration.yaml",
    update=True,
)


samples_PNetReg15_dict = {
    "2022_preEE": "QCD_PT-15to7000_PNetReg15_JMENano_Summer22",
    "2022_postEE": "QCD_PT-15to7000_PNetReg15_JMENano_Summer22EE",
    "2023_preBPix": "QCD_PT-15to7000_PNetReg15_JMENano_Summer23",
    "2023_postBPix": "QCD_PT-15to7000_PNetReg15_JMENano_Summer23BPix",
}


common_cats = {
    "baseline": [passthrough],
}

cfg = Configurator(
    parameters=parameters,
    datasets={
        "jsons": [
            # f"{localdir}/datasets/QCD.json",
            f"{localdir}/datasets/QCD_PNetReg15.json",
        ],
        "filter": {
            "samples": [
                (
                    samples_PNetReg15_dict[year]
                )
            ],
            "samples_exclude": [],
            "year": [year],
        },
        "subsamples": {},
    },
    workflow=METProcessor,
    workflow_options={
        "donotscale_sumgenweights": True,
        "only_physisical_jet":False,
    },
    skim=[],
    preselections=[PV_presel],
    categories={
                **common_cats,
    },
    weights={
        "common": {
            "inclusive": [
                # "genWeight",
                # "lumi",
                # "XS",
            ],
            "bycategory": {},
        },
        "bysample": {},
    },
    variations={
        "weights": {
            "common": {
                "inclusive": [],
                "bycategory": {},
            },
            "bysample": {},
        }
    },
    variables={
        **met_hists("PuppiMET"),
        **met_hists("PuppiMETPNet"),
        **met_hists("PuppiMETPNetPlusNeutrino"),
        **met_hists("GenMET"),
        **met_hists("GenMETPlusNeutrino"),
    },
    columns={
        "common": {
            "inclusive": [
                ColOut("PuppiMET", [ "pt", "phi"]),
                ColOut("PuppiMETPNet", [ "pt", "phi"]),
                ColOut("PuppiMETPNetPlusNeutrino", [ "pt", "phi"]),
                ColOut("GenMET", [ "pt", "phi"]),
                ColOut("GenMETPlusNeutrino", [ "pt", "phi"]),
            ]
        },
        "bysample": {
        },
    },
)
