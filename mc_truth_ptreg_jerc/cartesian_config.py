from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.parameters.cuts import passthrough
from pocket_coffea.lib.columns_manager import ColOut
from pocket_coffea.lib.categorization import CartesianSelection, MultiCut

from pocket_coffea.parameters import defaults
import os

from mc_truth_ptreg_jerc.workflow import *
from mc_truth_ptreg_jerc.cuts import *
from mc_truth_ptreg_jerc.custom_functions import get_closure_function_information
from mc_truth_ptreg_jerc.params.l2relative_txt import get_corr_file


from mc_truth_ptreg_jerc.custom_cut_functions import *
from mc_truth_ptreg_jerc.params.binning import *
from mc_truth_ptreg_jerc.variables_def import get_variables_dict


localdir = os.path.dirname(os.path.abspath(__file__))

# Loading default parameters

default_parameters = defaults.get_default_parameters()
defaults.register_configuration_dir("config_dir", localdir + "/params")


year = os.environ.get("YEAR", "2022_preEE")
# adding object preselection
parameters = defaults.merge_parameters_from_files(
    default_parameters,
    f"{localdir}/params/object_preselection{'_extendedPT' if (int(os.environ.get('EXTENDED_PT_BINS', 0)) == 1) else ''}.yaml",
    f"{localdir}/params/jets_calibration.yaml",
    update=True,
)

if int(os.environ.get("UPART", 0)) == 0:
    eta_bins = eta_bins
else:
    eta_bins = eta_bins_upart



# Loading the L2 Relative (MC truth) corrections
param_path = localdir + "/params"
(
    corr_files_pnetreg,
    corr_files_upartreg,
    corr_files_pnetreg_neutrino,
    corr_files_upartreg_neutrino,
    corr_files,
) = get_corr_file(param_path)

if int(os.environ.get("CLOSURE", 0)) == 1:

    print(f"Performing closure test with {corr_files_pnetreg[year]}")
    mc_truth_corr_pnetreg = get_closure_function_information(corr_files_pnetreg[year])

    print(f"Performing closure test with {corr_files_upartreg[year]}")
    mc_truth_corr_upartreg = get_closure_function_information(corr_files_upartreg[year])

    print(f"Performing closure test with {corr_files_pnetreg_neutrino[year]}")
    mc_truth_corr_pnetreg_neutrino = get_closure_function_information(
        corr_files_pnetreg_neutrino[year]
    )

    print(f"Performing closure test with {corr_files_upartreg_neutrino[year]}")
    mc_truth_corr_upartreg_neutrino = get_closure_function_information(
        corr_files_upartreg_neutrino[year]
    )

else:
    mc_truth_corr_pnetreg = None
    mc_truth_corr_pnetreg_neutrino = None
    mc_truth_corr_upartreg = None
    mc_truth_corr_upartreg_neutrino = None

print(f"Reapplying corrections {corr_files[year]}")
mc_truth_corr = get_closure_function_information(corr_files[year])


cuts_eta = []
cuts_names_eta = []
cuts_eta_neutrino = []
cuts_names_eta_neutrino = []
cuts_reco_eta = []
cuts_names_reco_eta = []

if int(os.environ.get("NEUTRINO", 1)) == 0:
    print("RECO JET ETA CUTS NEUTRINO==0")
    for i in range(len(eta_bins) - 1):
        eta_low, eta_high = eta_bins[i], eta_bins[i + 1]
        cuts_reco_eta.append(get_reco_etabin(eta_low, eta_high))
        cuts_names_reco_eta.append(f"MatchedJets_reco_eta{eta_low}to{eta_high}")
elif int(os.environ.get("NEUTRINO", 0)) == 1:
    print("RECO JET ETA CUTS NEUTRINO==1")
    for i in range(len(eta_bins) - 1):
        eta_low, eta_high = eta_bins[i], eta_bins[i + 1]
        cuts_reco_eta.append(get_reco_neutrino_etabin(eta_low, eta_high))
        cuts_names_reco_eta.append(f"MatchedJetsNeutrino_reco_eta{eta_low}to{eta_high}")
elif int(os.environ.get("ABS_ETA_INCLUSIVE", 0)) == 1:
    print("RUNNING ABS_ETA_INCLUSIVE ETA BINS")
    for i in range(len(eta_bins) - 1):
        eta_low, eta_high = eta_bins[i], eta_bins[i + 1]
        cuts_reco_eta.append(get_reco_neutrino_abs_etabin(eta_low, eta_high))
        cuts_names_reco_eta.append(
            f"MatchedJetsNeutrino_reco_abseta{eta_low}to{eta_high}"
        )
else:
    print("RECO JET ETA CUTS")
    for i in range(len(eta_bins) - 1):
        eta_low, eta_high = eta_bins[i], eta_bins[i + 1]
        cuts_reco_eta.append(get_reco_neutrino_etabin(eta_low, eta_high))
        cuts_names_reco_eta.append(f"MatchedJetsNeutrino_reco_eta{eta_low}to{eta_high}")


variables_dict = get_variables_dict(
    cuts_names_eta, cuts_names_reco_eta, cuts_names_eta_neutrino
)

# samples_dict = {
# "2022_preEE": "QCD_PT-15to7000_JMENano_Summer22",
# "2022_postEE": "QCD_PT-15to7000_JMENano_Summer22EE",
# "2023_preBPix": "QCD_PT-15to7000_JMENano_Summer23",
# "2023_postBPix": "QCD_PT-15to7000_JMENano_Summer23BPix",
# }

samples_dict = {
    "2016_PreVFP": "QCD_PT-15to7000_JMENano_Summer20UL16APV",
    "2016_PostVFP": "QCD_PT-15to7000_JMENano_Summer20UL16",
    "2017": "QCD_PT-15to7000_JMENano_Summer20UL17",
    "2018": "QCD_PT-15to7000_JMENano_Summer20UL18",
    "2022_preEE": "QCD_PT-15to7000_JMENanov15_Summer22",
    "2022_postEE": "QCD_PT-15to7000_JMENanov15_Summer22EE",
    "2023_preBPix": "QCD_PT-15to7000_JMENanov15_Summer23",
    "2023_postBPix": "QCD_PT-15to7000_JMENanov15_Summer23BPix",
    "2024": "QCD_PT-15to7000_JMENano_Summer24",
    "2025": "QCD_PT-15to7000_JMENano_Winter25",
}
samples_PNetReg15_dict = {
    "2022_preEE": "QCD_PT-15to7000_PNetReg15_JMENano_Summer22",
    "2022_postEE": "QCD_PT-15to7000_PNetReg15_JMENano_Summer22EE",
    "2023_preBPix": "QCD_PT-15to7000_PNetReg15_JMENano_Summer23",
    "2023_postBPix": "QCD_PT-15to7000_PNetReg15_JMENano_Summer23BPix",
}
# No UparTReg15 used for now but may be added later if necessary

multicuts = [
    MultiCut(
        name="eta",
        cuts=cuts_eta + cuts_eta_neutrino + cuts_reco_eta,
        cuts_names=cuts_names_eta + cuts_names_eta_neutrino + cuts_names_reco_eta,
    ),
]

common_cats = {
    # "baseline": [passthrough],
    "PV_presel": [PV_presel],
}

cfg = Configurator(
    parameters=parameters,
    datasets={
        "jsons": [
            # f"{localdir}/datasets/QCD.json",
            f"{localdir}/datasets/QCD_redirector.json",
            f"{localdir}/datasets/QCD_PNetReg15.json",
        ],
        "filter": {
            "samples": [
                (
                    samples_PNetReg15_dict[year]
                    if (
                        int(os.environ.get("PNETREG15", 0)) == 1
                        or int(os.environ.get("SPLITPNETREG15", 0)) == 1
                    )
                    else samples_dict[year]
                )
            ],
            "samples_exclude": [],
            "year": [year],
        },
        "subsamples": {},
    },
    workflow=QCDBaseProcessor,
    workflow_options={
        "donotscale_sumgenweights": True,
        "mc_truth_corr_pnetreg": mc_truth_corr_pnetreg,
        "mc_truth_corr_pnetreg_neutrino": mc_truth_corr_pnetreg_neutrino,
        "mc_truth_corr_upartreg": mc_truth_corr_upartreg,
        "mc_truth_corr_upartreg_neutrino": mc_truth_corr_upartreg_neutrino,
        "mc_truth_corr": mc_truth_corr,
        "DeltaR_matching": 0.2,
        "SetRegResponseToZero": True,
        "GenJetPtCut": (
            15
            if (
                int(os.environ.get("PNETREG15", 0)) == 1
                or int(os.environ.get("SPLITPNETREG15", 0)) == 1
            )
            else (0 if (int(os.environ.get("EXTENDED_PT_BINS", 0)) == 1) else 50)
        ),
        "pnet": False,
        "upart": False,
    },
    calibrators=[],
    skim=[],
    preselections=[],
    categories=CartesianSelection(multicuts=multicuts, common_cats=common_cats),
    # categories={
    #             **common_cats,
    #     # **eta_cuts,
    # },
    weights={
        "common": {
            "inclusive": [],
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
    variables=variables_dict,
    columns={},
)
