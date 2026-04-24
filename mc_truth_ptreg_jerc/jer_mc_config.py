import os
import hist

from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.parameters.cuts import passthrough
from pocket_coffea.lib.columns_manager import ColOut
from pocket_coffea.lib.hist_manager import HistConf, Axis
from pocket_coffea.lib.categorization import CartesianSelection, MultiCut

from pocket_coffea.parameters import defaults

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

DUMP_COLUMNS_AS_ARRAYS_PER_CHUNK = False

year = os.environ.get("YEAR", "2022_preEE")
print("YEAR", year)
add_str = ""
output_chunks_name = f"root://t3dcachedb03.psi.ch:1094//pnfs/psi.ch/cms/trivcat/store/user/mmalucch/out_jer_mc_ptreg/out_jer_closure_pnet_upart_{year}{add_str}/parquet_files"
if DUMP_COLUMNS_AS_ARRAYS_PER_CHUNK:
    print("Output chunks path:", output_chunks_name)


# adding object preselection
parameters = defaults.merge_parameters_from_files(
    default_parameters,
    f"{localdir}/params/object_preselection{'_extendedPT' if (int(os.environ.get('EXTENDED_PT_BINS', 0)) == 1) else ''}.yaml",
    f"{localdir}/params/jets_calibration.yaml",
    update=True,
)

# Loading the L2 Relative (MC truth) corrections
param_path = localdir + "/params"
(
    corr_files_pnetreg,
    corr_files_upartreg,
    corr_files_pnetreg_neutrino,
    corr_files_upartreg_neutrino,
    corr_files,
) = get_corr_file(param_path)

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

print(f"Reapplying corrections {corr_files[year]}")
mc_truth_corr = get_closure_function_information(corr_files[year])


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


cfg = Configurator(
    parameters=parameters,
    datasets={
        "jsons": [
            f"{localdir}/datasets/QCD.json",
            # f"{localdir}/datasets/QCD_redirector.json",
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
        "pnet": True,
        "upart": True,
        "dump_columns_as_arrays_per_chunk": (
            output_chunks_name if DUMP_COLUMNS_AS_ARRAYS_PER_CHUNK else None
        ),
    },
    calibrators=[],
    skim=[],
    preselections=[PV_presel],
    categories={
        "baseline": [passthrough],
    },
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
    variables={
        "rho_pu": HistConf(
            [
                Axis(
                    coll="Pileup",
                    field="nPU",
                    bins=PU_bins,
                    label="Pileup_nPU",
                ),
                Axis(
                    coll="Rho",
                    field="fixedGridRhoFastjetAll",
                    bins=rho_bins,
                    label="Rho_fixedGridRhoFastjetAll",
                ),
            ],
            storage="mean",
        )
    },
    columns={
        "common": {
            "inclusive": (
                (
                    [
                        ColOut("Rho", ["fixedGridRhoFastjetAll"]),
                        ColOut("Pileup", ["nPU"]),
                        ColOut(
                            "MatchedJets",
                            [
                                "pt",
                                "eta",
                                "ResponseJEC",
                                "JetPtJEC",
                                "ResponseRaw",
                                "JetPtRaw",
                                "ResponsePNetReg",
                                "JetPtPNetReg",
                                "ResponseUparTReg",
                                "JetPtUparTReg",
                            ],
                            flatten=False if DUMP_COLUMNS_AS_ARRAYS_PER_CHUNK else True,
                        ),
                        ColOut(
                            "MatchedJetsNeutrino",
                            [
                                "pt",
                                "eta",
                                "ResponseJEC",
                                "JetPtJEC",
                                "ResponseRaw",
                                "JetPtRaw",
                                "ResponsePNetRegNeutrino",
                                "JetPtPNetRegNeutrino",
                                "ResponseUparTRegNeutrino",
                                "JetPtUparTRegNeutrino",
                            ],
                            flatten=False if DUMP_COLUMNS_AS_ARRAYS_PER_CHUNK else True,
                        ),
                    ]
                )
            ),
        },
    },
)
