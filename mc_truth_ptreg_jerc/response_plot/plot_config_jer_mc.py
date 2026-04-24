import numpy as np

TEST = False

MIXED_MODE = True

REBIN_FOR_PLOTTING = False

pu_bins = np.array(
    [
        0,
        10,
        20,
        30,
        40,
        50,
        60,
        70,
        80,
        90,
        100,
    ]
)

jet_eta_bins = np.array(
    [
        -5.191,
        -3.839,
        -3.489,
        -3.139,
        -2.964,
        -2.853,
        -2.65,
        -2.5,
        -2.322,
        -2.172,
        -2.043,
        -1.93,
        -1.74,
        -1.566,
        -1.305,
        -1.044,
        -0.783,
        -0.522,
        -0.261,
        0.0,
        0.261,
        0.522,
        0.783,
        1.044,
        1.305,
        1.566,
        1.74,
        1.93,
        2.043,
        2.172,
        2.322,
        2.5,
        2.65,
        2.853,
        2.964,
        3.139,
        3.489,
        3.839,
        5.191,
    ]
)

jet_pt_bins = np.array(
    [
        8,
        10,
        12,
        15,
        17,
        20,
        23,
        27,
        30,
        35,
        40,
        45,
        57,
        72,
        90,
        120,
        150,
        200,
        300,
        400,
        550,
        750,
        1000,
        1500,
        2000,
        2500,
        3000,
        3500,
        4000,
        4500,
        5000,
    ]
)

if TEST:
    pu_bins = np.array([0, 50, 100])
    jet_eta_bins = np.array([-5.191, 0.0, 5.191])
    jet_pt_bins = np.array(
        [8, 50, 120, 200, 300, 400, 550, 750, 1000, 1500, 2000, 2500, 3000, 4000, 5000]
    )

BIN_VARIABLES = {
    "Pileup_nPU": {
        # if True, expand the variable to have it at jet-level
        "event_var": True,
        "bin_edges": pu_bins,
        "label": r"$\mu$",
        "name_plot": "nPU",
        "resolution_x_variable": False,
        "txt_map_to": "Rho_fixedGridRhoFastjetAll",
    },
    "MatchedJets_eta": {
        "event_var": False,
        "bin_edges": jet_eta_bins,
        "label": r"$\eta$",
        "name_plot": "jet_eta",
        "resolution_x_variable": False,
        "txt_name": "JetEta",
    },
    "MatchedJets_pt": {
        "event_var": False,
        "bin_edges": jet_pt_bins,
        "label": r"$p_{\mathrm{T}}^{ptcl}$ [GeV]",
        "name_plot": "jet_gen_pt",
        "resolution_x_variable": True,
        "txt_name": "JetPt",
    },
}

BIN_VARIABLES_NEUTRINO = {
    "Pileup_nPU": {
        # if True, expand the variable to have it at jet-level
        "event_var": True,
        "bin_edges": pu_bins,
        "label": r"$\mu$",
        "name_plot": "nPU",
        "resolution_x_variable": False,
        "txt_map_to": "Rho_fixedGridRhoFastjetAll",
    },
    "MatchedJetsNeutrino_eta": {
        "event_var": False,
        "bin_edges": jet_eta_bins,
        "label": r"$\eta$",
        "name_plot": "jet_eta",
        "resolution_x_variable": False,
        "txt_name": "JetEta",
    },
    "MatchedJetsNeutrino_pt": {
        "event_var": False,
        "bin_edges": jet_pt_bins,
        "label": r"$p_{\mathrm{T}}^{ptcl}$ [GeV]",
        "name_plot": "jet_gen_pt",
        "resolution_x_variable": True,
        "txt_name": "JetPt",
    },
}


MAPPING_VARIABLES = {
    "Rho_fixedGridRhoFastjetAll": {
        "event_var": True,
        "label": r"$\rho$",
        "name_plot": "rho",
        "N_bins": 1000,
        "bin_vars": ["Pileup_nPU"],
        "linear_fit": True,
        "txt_name": "Rho",
        "rebin_for_plotting": True,
    },
    "MatchedJets_JetPtJEC": {
        "label": r"$p_{\mathrm{T}}^{\mathrm{reco}}$ [GeV]",
        "name_plot": "jet_reco_pt",
        "N_bins": 1000,
        "bin_limits": (0, 6000),
        "bin_vars": ["Pileup_nPU", "MatchedJets_eta", "MatchedJets_pt"],
        "legend_name": "JEC",
        "rebin_for_plotting": True,
    },
    "MatchedJets_JetPtRaw": {
        "label": r"$p_{\mathrm{T}}^{\mathrm{reco}}$ [GeV]",
        "name_plot": "jet_reco_pt",
        "N_bins": 1000,
        "bin_limits": (0, 6000),
        "bin_vars": ["Pileup_nPU", "MatchedJets_eta", "MatchedJets_pt"],
        "legend_name": "Raw",
        "rebin_for_plotting": True,
    },
    "MatchedJets_JetPtPNetReg": {
        "label": r"$p_{\mathrm{T}}^{\mathrm{reco}}$ [GeV]",
        "name_plot": "jet_reco_pt",
        "N_bins": 1000,
        "bin_limits": (0, 6000),
        "bin_vars": ["Pileup_nPU", "MatchedJets_eta", "MatchedJets_pt"],
        "legend_name": "PNetReg",
        "rebin_for_plotting": True,
    },
    "MatchedJets_JetPtUparTReg": {
        "label": r"$p_{\mathrm{T}}^{\mathrm{reco}}$ [GeV]",
        "name_plot": "jet_reco_pt",
        "N_bins": 1000,
        "bin_limits": (0, 6000),
        "bin_vars": ["Pileup_nPU", "MatchedJets_eta", "MatchedJets_pt"],
        "legend_name": "UparTReg",
        "rebin_for_plotting": True,
    },
    "MatchedJetsNeutrino_JetPtJEC": {
        "label": r"$p_{\mathrm{T}}^{\mathrm{reco}}$ [GeV]",
        "name_plot": "jet_reco_pt_neutrino",
        "N_bins": 1000,
        "bin_limits": (0, 6000),
        "bin_vars": ["Pileup_nPU", "MatchedJetsNeutrino_eta", "MatchedJetsNeutrino_pt"],
        "legend_name": "JEC",
        "rebin_for_plotting": True,
    },
    "MatchedJetsNeutrino_JetPtRaw": {
        "label": r"$p_{\mathrm{T}}^{\mathrm{reco}}$ [GeV]",
        "name_plot": "jet_reco_pt_neutrino",
        "N_bins": 1000,
        "bin_limits": (0, 6000),
        "bin_vars": ["Pileup_nPU", "MatchedJetsNeutrino_eta", "MatchedJetsNeutrino_pt"],
        "legend_name": "Raw",
        "rebin_for_plotting": True,
    },
    "MatchedJetsNeutrino_JetPtPNetRegNeutrino": {
        "label": r"$p_{\mathrm{T}}^{\mathrm{reco}}$ [GeV]",
        "name_plot": "jet_reco_pt_neutrino",
        "N_bins": 1000,
        "bin_limits": (0, 6000),
        "bin_vars": ["Pileup_nPU", "MatchedJetsNeutrino_eta", "MatchedJetsNeutrino_pt"],
        "legend_name": "PNetRegNeutrino",
        "rebin_for_plotting": True,
    },
    "MatchedJetsNeutrino_JetPtUparTRegNeutrino": {
        "label": r"$p_{\mathrm{T}}^{\mathrm{reco}}$ [GeV]",
        "name_plot": "jet_reco_pt_neutrino",
        "N_bins": 1000,
        "bin_limits": (0, 6000),
        "bin_vars": ["Pileup_nPU", "MatchedJetsNeutrino_eta", "MatchedJetsNeutrino_pt"],
        "legend_name": "UparTRegNeutrino",
        "rebin_for_plotting": True,
    },
}

RESPONSE_VARIABLES = {
    "MatchedJets_ResponseJEC": {
        "label": r"$p_{\mathrm{T}}^{\mathrm{reco}}$/$p_{\mathrm{T}}^{ptcl}$",
        "name_plot": "response_jec",
        "N_bins": 60,
        "bin_limits": (0, 2.0),
        "bin_vars": ["Pileup_nPU", "MatchedJets_eta", "MatchedJets_pt"],
        "legend_name": "JEC",
        "map_x_variable": "MatchedJets_JetPtJEC",
        "is_reference": True,
        "txt_jet_name": "AK4PFPuppi",
        "rebin_for_plotting": False,
        "rebin_for_plotting": False,
    },
    "MatchedJets_ResponseRaw": {
        "label": r"$p_{\mathrm{T}}^{\mathrm{reco}}$/$p_{\mathrm{T}}^{ptcl}$",
        "name_plot": "response_raw",
        "N_bins": 60,
        "bin_limits": (0, 2.0),
        "bin_vars": ["Pileup_nPU", "MatchedJets_eta", "MatchedJets_pt"],
        "legend_name": "Raw",
        "map_x_variable": "MatchedJets_JetPtRaw",
        "is_reference": False,
        "txt_jet_name": "AK4PFPuppiRaw",
        "rebin_for_plotting": False,
    },
    "MatchedJets_ResponsePNetReg": {
        "label": r"$p_{\mathrm{T}}^{\mathrm{reco}}$/$p_{\mathrm{T}}^{ptcl}$",
        "name_plot": "response_pnetreg",
        "N_bins": 60,
        "bin_limits": (0, 2.0),
        "bin_vars": ["Pileup_nPU", "MatchedJets_eta", "MatchedJets_pt"],
        "legend_name": "PNetReg",
        "map_x_variable": "MatchedJets_JetPtPNetReg",
        "is_reference": False,
        "txt_jet_name": "AK4PFPNet",
        "rebin_for_plotting": False,
    },
    "MatchedJets_ResponseUparTReg": {
        "label": r"$p_{\mathrm{T}}^{\mathrm{reco}}$/$p_{\mathrm{T}}^{ptcl}$",
        "name_plot": "response_upartreg",
        "N_bins": 60,
        "bin_limits": (0, 2.0),
        "bin_vars": ["Pileup_nPU", "MatchedJets_eta", "MatchedJets_pt"],
        "legend_name": "UParTReg",
        "map_x_variable": "MatchedJets_JetPtUparTReg",
        "is_reference": False,
        "txt_jet_name": "AK4PFUParT",
        "rebin_for_plotting": False,
    },
}

RESPONSE_VARIABLES_NEUTRINO = {
    "MatchedJetsNeutrino_ResponseJEC": {
        "label": r"$p_{\mathrm{T}}^{\mathrm{reco}}$/$p_{\mathrm{T}}^{ptcl}$",
        "name_plot": "response_jec",
        "N_bins": 60,
        "bin_limits": (0, 2.0),
        "bin_vars": ["Pileup_nPU", "MatchedJetsNeutrino_eta", "MatchedJetsNeutrino_pt"],
        "legend_name": "JEC",
        "map_x_variable": "MatchedJetsNeutrino_JetPtJEC",
        "is_reference": True,
        "txt_jet_name": "AK4PFPuppiNeutrino",
        "rebin_for_plotting": False,
    },
    "MatchedJetsNeutrino_ResponseRaw": {
        "label": r"$p_{\mathrm{T}}^{\mathrm{reco}}$/$p_{\mathrm{T}}^{ptcl}$",
        "name_plot": "response_raw",
        "N_bins": 60,
        "bin_limits": (0, 2.0),
        "bin_vars": ["Pileup_nPU", "MatchedJetsNeutrino_eta", "MatchedJetsNeutrino_pt"],
        "legend_name": "Raw",
        "map_x_variable": "MatchedJetsNeutrino_JetPtRaw",
        "is_reference": False,
        "txt_jet_name": "AK4PFPuppiNeutrinoRaw",
        "rebin_for_plotting": False,
    },
    "MatchedJetsNeutrino_ResponsePNetRegNeutrino": {
        "label": r"$p_{\mathrm{T}}^{\mathrm{reco}}$/$p_{\mathrm{T}}^{ptcl}$",
        "name_plot": "response_pnetreg_neutrino",
        "N_bins": 60,
        "bin_limits": (0, 2.0),
        "bin_vars": ["Pileup_nPU", "MatchedJetsNeutrino_eta", "MatchedJetsNeutrino_pt"],
        "legend_name": "PNetRegNeutrino",
        "map_x_variable": "MatchedJetsNeutrino_JetPtPNetRegNeutrino",
        "is_reference": False,
        "txt_jet_name": "AK4PFPNetPlusNeutrino",
        "rebin_for_plotting": False,
    },
    "MatchedJetsNeutrino_ResponseUparTRegNeutrino": {
        "label": r"$p_{\mathrm{T}}^{\mathrm{reco}}$/$p_{\mathrm{T}}^{ptcl}$",
        "name_plot": "response_upartreg_neutrino",
        "N_bins": 60,
        "bin_limits": (0, 2.0),
        "bin_vars": ["Pileup_nPU", "MatchedJetsNeutrino_eta", "MatchedJetsNeutrino_pt"],
        "legend_name": "UParTRegNeutrino",
        "map_x_variable": "MatchedJetsNeutrino_JetPtUparTRegNeutrino",
        "is_reference": False,
        "txt_jet_name": "AK4PFUParTPlusNeutrino",
        "rebin_for_plotting": False,
    },
}


RESPONSE_VARIABLES_MIXED = {
    "MatchedJets_ResponseJEC": {
        "label": r"$p_{\mathrm{T}}^{\mathrm{reco}}$/$p_{\mathrm{T}}^{ptcl}$",
        "name_plot": "response_jec",
        "N_bins": 60,
        "bin_limits": (0, 2.0),
        "bin_vars": ["Pileup_nPU", "MatchedJets_eta", "MatchedJets_pt"],
        "legend_name": "JEC",
        "map_x_variable": "MatchedJets_JetPtJEC",
        "is_reference": True,
        "txt_jet_name": "AK4PFPuppi",
        "rebin_for_plotting": False,
    },
    "MatchedJets_ResponseRaw": {
        "label": r"$p_{\mathrm{T}}^{\mathrm{reco}}$/$p_{\mathrm{T}}^{ptcl}$",
        "name_plot": "response_raw",
        "N_bins": 60,
        "bin_limits": (0, 2.0),
        "bin_vars": ["Pileup_nPU", "MatchedJets_eta", "MatchedJets_pt"],
        "legend_name": "Raw",
        "map_x_variable": "MatchedJets_JetPtRaw",
        "is_reference": False,
        "txt_jet_name": "AK4PFPuppiRaw",
        "rebin_for_plotting": False,
    },
    "MatchedJets_ResponsePNetReg": {
        "label": r"$p_{\mathrm{T}}^{\mathrm{reco}}$/$p_{\mathrm{T}}^{ptcl}$",
        "name_plot": "response_pnetreg",
        "N_bins": 60,
        "bin_limits": (0, 2.0),
        "bin_vars": ["Pileup_nPU", "MatchedJets_eta", "MatchedJets_pt"],
        "legend_name": "PNetReg",
        "map_x_variable": "MatchedJets_JetPtPNetReg",
        "is_reference": False,
        "txt_jet_name": "AK4PFPNet",
        "rebin_for_plotting": False,
    },
    "MatchedJets_ResponseUparTReg": {
        "label": r"$p_{\mathrm{T}}^{\mathrm{reco}}$/$p_{\mathrm{T}}^{ptcl}$",
        "name_plot": "response_upartreg",
        "N_bins": 60,
        "bin_limits": (0, 2.0),
        "bin_vars": ["Pileup_nPU", "MatchedJets_eta", "MatchedJets_pt"],
        "legend_name": "UParTReg",
        "map_x_variable": "MatchedJets_JetPtUparTReg",
        "is_reference": False,
        "txt_jet_name": "AK4PFUParT",
        "rebin_for_plotting": False,
    },
    "MatchedJetsNeutrino_ResponsePNetRegNeutrino": {
        "label": r"$p_{\mathrm{T}}^{\mathrm{reco}}$/$p_{\mathrm{T}}^{ptcl}$",
        "name_plot": "response_pnetreg_neutrino",
        "N_bins": 60,
        "bin_limits": (0, 2.0),
        "bin_vars": ["Pileup_nPU", "MatchedJetsNeutrino_eta", "MatchedJetsNeutrino_pt"],
        "legend_name": "PNetRegNeutrino",
        "map_x_variable": "MatchedJetsNeutrino_JetPtPNetRegNeutrino",
        "is_reference": False,
        "txt_jet_name": "AK4PFPNetPlusNeutrino",
        "rebin_for_plotting": False,
    },
    "MatchedJetsNeutrino_ResponseUparTRegNeutrino": {
        "label": r"$p_{\mathrm{T}}^{\mathrm{reco}}$/$p_{\mathrm{T}}^{ptcl}$",
        "name_plot": "response_upartreg_neutrino",
        "N_bins": 60,
        "bin_limits": (0, 2.0),
        "bin_vars": ["Pileup_nPU", "MatchedJetsNeutrino_eta", "MatchedJetsNeutrino_pt"],
        "legend_name": "UParTRegNeutrino",
        "map_x_variable": "MatchedJetsNeutrino_JetPtUparTRegNeutrino",
        "is_reference": False,
        "txt_jet_name": "AK4PFUParTPlusNeutrino",
        "rebin_for_plotting": False,
    },
}

# Merged dictionaries combining regular and neutrino versions
BIN_VARIABLES_MIXED = {
    **BIN_VARIABLES,
    **BIN_VARIABLES_NEUTRINO,
}

PLOT_SETTINGS_DICT = {
    "JEC": {"color": "darkorange", "fmt": "o"},
    "Raw": {"color": "green", "fmt": "s"},
    "PNetReg": {"color": "darkred", "fmt": "<"},
    "UparTReg": {"color": "tomato", "fmt": "^"},
    "PNetRegNeutrino": {"color": "darkblue", "fmt": ">"},
    "UparTRegNeutrino": {"color": "cornflowerblue", "fmt": "v"},
}

Y_LIM_RESOLUTION = (0, 0.35)

YEAR_MAP = {
    "2022_preEE": "Summer22",
    "2022_postEE": "Summer22EE",
    "2023_preBPix": "Summer23",
    "2023_postBPix": "Summer23BPix",
    "2024": "Winter24",
}
