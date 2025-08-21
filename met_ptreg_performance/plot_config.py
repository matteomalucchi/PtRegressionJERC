import numpy as np
import mplhep as hep

# var_dict = {
#     r"MET $p_{\mathrm{T}}$ [GeV]": [
#         "PuppiMET_pt",
#         "PuppiMETPNet_pt",
#         "PuppiMETPNetPlusNeutrino_pt",
#     ],
#     r"MET $\phi$": ["PuppiMET_phi", "PuppiMETPNet_phi", "PuppiMETPNetPlusNeutrino_phi"],
#     r"MET_MinusMuons $p_{\mathrm{T}}$ [GeV]": [
#         "PuppiMET_MuonGood_pt",
#         "PuppiMETPNet_MuonGood_pt",
#         "PuppiMETPNetPlusNeutrino_MuonGood_pt",
#     ],
#     r"MET_MinusMuons $\phi$": [
#         "PuppiMET_MuonGood_phi",
#         "PuppiMETPNet_MuonGood_phi",
#         "PuppiMETPNetPlusNeutrino_MuonGood_phi",
#     ],
#     r"$-u_{||}/q_{T}$": [
#         "PuppiMET_MuonGood_response",
#         "PuppiMETPNet_MuonGood_response",
#         "PuppiMETPNetPlusNeutrino_MuonGood_response",
#     ],
#     r"$u_{||}+q_{T}$ [GeV]": [
#         "PuppiMET_MuonGood_u_paral_predict",
#         "PuppiMETPNet_MuonGood_u_paral_predict",
#         "PuppiMETPNetPlusNeutrino_MuonGood_u_paral_predict",
#     ],
#     r"$u_{\perp}$ [GeV]": [
#         "PuppiMET_MuonGood_u_perp_predict",
#         "PuppiMETPNet_MuonGood_u_perp_predict",
#         "PuppiMETPNetPlusNeutrino_MuonGood_u_perp_predict",
#     ],
#     r"MET $p_{\mathrm{T}}$ [GeV]": ["PuppiMET_pt", "PuppiMET_MuonGood_pt"],
#     r"METPNet $p_{\mathrm{T}}$ [GeV]": ["PuppiMETPNet_pt", "PuppiMETPNet_MuonGood_pt"],
#     r"METPNetPlusNeutrino $p_{\mathrm{T}}$ [GeV]": [
#         "PuppiMETPNetPlusNeutrino_pt",
#         "PuppiMETPNetPlusNeutrino_MuonGood_pt",
#     ],
# }
# ranges = {
#     r"MET $p_{\mathrm{T}}$ [GeV]": (0, 200),
#     r"MET $\phi$": (-3.14, 3.14),
#     r"MET_MinusMuons $p_{\mathrm{T}}$ [GeV]": (0, 200),
#     r"MET_MinusMuons $\phi$": (-3.14, 3.14),
#     r"$-u_{||}/q_{T}$": (-1, 2),
#     r"$u_{||}+q_{T}$ [GeV]": (-200, 200),
#     r"$u_{\perp}$ [GeV]": (-200, 200),
#     r"MET $p_{\mathrm{T}}$ [GeV]": (0, 200),
#     r"METPNet $p_{\mathrm{T}}$ [GeV]": (0, 200),
#     r"METPNetPlusNeutrino $p_{\mathrm{T}}$ [GeV]": (0, 200),
# }
# log_dict = {
#     r"MET $p_{\mathrm{T}}$ [GeV]": True,
#     r"MET $\phi$": False,
#     r"MET_MinusMuons $p_{\mathrm{T}}$ [GeV]": True,
#     r"MET_MinusMuons $\phi$": False,
#     r"$-u_{||}/q_{T}$": True,
#     r"$u_{||}+q_{T}$ [GeV]": True,
#     r"$u_{\perp}$ [GeV]": True,
#     r"MET $p_{\mathrm{T}}$ [GeV]": True,
#     r"METPNet $p_{\mathrm{T}}$ [GeV]": True,
#     r"METPNetPlusNeutrino $p_{\mathrm{T}}$ [GeV]": True,
# }
# color_list = [
#     "black",
#     "red",
#     "green",
#     "orange",
#     "purple",
#     "brown",
#     "pink",
#     "gray",
#     "olive",
#     "cyan",
# ]


met_list = [
    "RawPuppiMET",
    "RawPuppiMET-Type1",
    "RawPuppiMET-Type1JEC",
    "RawPuppiMET-Type1PNet",
    "RawPuppiMET-Type1PNetPlusNeutrino",
    "PuppiMET",
    "PuppiMET-Type1",
    "PuppiMET-Type1JEC",
    "PuppiMET-Type1PNet",
    "PuppiMET-Type1PNetPlusNeutrino",
]

hadronic_recoil_dict={
    f"{met}_pt": {
        "plot_name": met+r" $p_{\mathrm{T}}$ [GeV]",
        "variables": [f"{met}_pt", f"u{met}_pt"],
        "range": (0, 200),
        "log": True,
        "ratio_label": "u / MET",
    } for met in met_list
}

total_var_dict = {
    "MET_comparison_pt": {
        "plot_name": r"MET $p_{\mathrm{T}}$ [GeV]",
        "variables": [met + "_pt" for met in met_list],
        "range": (0, 200),
        "log": True,
        "ratio_label": "MET / RawPuppiMET",
    },
    "MET_comparison_phi": {
        "plot_name": r"MET $\phi$",
        "variables": [met + "_phi" for met in met_list],
        "range": (-3.14, 3.14),
        "log": False,
        "ratio_label": "MET / RawPuppiMET",
    },
    "hadronic_recoil_comparison_pt": {
        "plot_name": r"u $p_{\mathrm{T}}$ [GeV]",
        "variables": ["u"+met + "_pt" for met in met_list],
        "range": (0, 200),
        "log": True,
        "ratio_label": "MET / RawPuppiMET",
    },
    "hadronic_recoil_comparison_phi": {
        "plot_name": r"u $\phi$",
        "variables": ["u"+met + "_phi" for met in met_list],
        "range": (-3.14, 3.14),
        "log": False,
        "ratio_label": "MET / RawPuppiMET",
    },
    "response_comparison": {
        "plot_name": r"$-u_{||}/q_{T}$",
        "variables": ["u"+met + "_response" for met in met_list],
        "range": (-1, 2),
        "log": True,
        "ratio_label": "MET / RawPuppiMET",
    },
    "u_paral_predict_comparison": {
        "plot_name": r"$u_{||}+q_{T}$ [GeV]",
        "variables": ["u"+met + "_u_paral_predict" for met in met_list],
        "range": (-200, 200),
        "log": True,
        "ratio_label": "MET / RawPuppiMET",
    },
    "u_perp_predict_comparison": {
        "plot_name": r"$u_{\perp}$ [GeV]",
        "variables": ["u"+met + "_u_perp_predict" for met in met_list],
        "range": (-200, 200),
        "log": True,
        "ratio_label": "MET / RawPuppiMET",
    },
    # "RawPuppiMET_pt": {
    #     "plot_name": r"RawPuppiMET $p_{\mathrm{T}}$ [GeV]",
    #     "variables": ["RawPuppiMET_pt", "hadronic_recoil_RawPuppiMET_pt"],
    #     "range": (0, 200),
    #     "log": True,
    #     "ratio_label": "u / MET",
    # },
    # "RawPuppiMETType1_pt": {
    #     "plot_name": r"RawPuppiMETType1 $p_{\mathrm{T}}$ [GeV]",
    #     "variables": ["RawPuppiMETType1_pt", "hadronic_recoil_RawPuppiMETType1_pt"],
    #     "range": (0, 200),
    #     "log": True,
    #     "ratio_label": "u / MET",
    # },
    # "RawPuppiMET_PNet_pt": {
    #     "plot_name": r"RawPuppiMET PNet $p_{\mathrm{T}}$ [GeV]",
    #     "variables": ["RawPuppiMETPNet_pt", "hadronic_recoil_RawPuppiMETPNet_pt"],
    #     "range": (0, 200),
    #     "log": True,
    #     "ratio_label": "u / MET",
    # },
    # "RawPuppiMET_PNetPlusNeutrino_pt": {
    #     "plot_name": r"RawPuppiMET PNetPlusNeutrino $p_{\mathrm{T}}$ [GeV]",
    #     "variables": [
    #         "RawPuppiMETPNetPlusNeutrino_pt",
    #         "hadronic_recoil_RawPuppiMETPNetPlusNeutrino_pt",
    #     ],
    #     "range": (0, 200),
    #     "log": True,
    #     "ratio_label": "u / MET",
    # },
    # "PuppiMET_pt": {
    #     "plot_name": r"PuppiMET $p_{\mathrm{T}}$ [GeV]",
    #     "variables": ["PuppiMET_pt", "hadronic_recoil_PuppiMET_pt"],
    #     "range": (0, 200),
    #     "log": True,
    #     "ratio_label": "u / MET",
    # },
    # "PuppiMETType1_pt": {
    #     "plot_name": r"PuppiMETType1 $p_{\mathrm{T}}$ [GeV]",
    #     "variables": ["PuppiMETType1_pt", "hadronic_recoil_PuppiMETType1_pt"],
    #     "range": (0, 200),
    #     "log": True,
    #     "ratio_label": "u / MET",
    # },
    # "PuppiMET_PNet_pt": {
    #     "plot_name": r"PuppiMET PNet $p_{\mathrm{T}}$ [GeV]",
    #     "variables": ["PuppiMETPNet_pt", "hadronic_recoil_PuppiMETPNet_pt"],
    #     "range": (0, 200),
    #     "log": True,
    #     "ratio_label": "u / MET",
    # },
    # "PuppiMET_PNetPlusNeutrino_pt": {
    #     "plot_name": r"PuppiMET PNetPlusNeutrino $p_{\mathrm{T}}$ [GeV]",
    #     "variables": [
    #         "PuppiMETPNetPlusNeutrino_pt",
    #         "hadronic_recoil_PuppiMETPNetPlusNeutrino_pt",
    #     ],
    #     "range": (0, 200),
    #     "log": True,
    #     "ratio_label": "u / MET",
    # },
}
total_var_dict.update(hadronic_recoil_dict)

response_var_name_dict = {
    "R": r"$-u_{\parallel}/q_{T}$",
    "R_mean": r"-$<u_{\parallel}/q_{T}>$",
    "R_quantile_resolution": r"$q \sigma(u_{\parallel}/q_{T})$",
    "R_stddev_resolution": r"$std dev \sigma(u_{\parallel}/q_{T})$",
    "u_perp": r"$u_{\perp}$ [GeV]",
    "u_perp_scaled": r"$u_{\perp}$ / (-$<u_{\parallel}/q_{T}>$)",
    "u_perp_mean": r"$<u_{\perp}>$ [GeV]",
    "u_perp_scaled_mean": r"$<u_{\perp}> / (-<u_{\parallel}/q_{T}>)$",
    "u_perp_quantile_resolution": r"$q \sigma(u_{\perp})$ [GeV]",
    "u_perp_scaled_quantile_resolution": r"$q \sigma(u_{\perp}) / (-<u_{\parallel}/q_{T}>)$",
    "u_perp_stddev_resolution": r"$std dev \sigma({u_{\perp}})$ [GeV]",
    "u_perp_stddev_scaled_resolution": r"$std dev \sigma({u_{\perp}}) / (-<u_{\parallel}/q_{T}>)$",
    "u_paral": r"$u_{\parallel}+q_{T}$ [GeV]",
    "u_paral_scaled": r"$u_{\parallel}+q_{T}$ / (-$<u_{\parallel}/q_{T}>$)",
    "u_paral_mean": r"$<u_{\parallel}>$ [GeV]",
    "u_paral_scaled_mean": r"$<u_{\parallel}> / (-<u_{\parallel}/q_{T}>)$",
    "u_paral_quantile_resolution": r"$q \sigma(u_{\parallel})$ [GeV]",
    "u_paral_scaled_quantile_resolution": r"$q \sigma(u_{\parallel}) / (-<u_{\parallel}/q_{T}>)$",
    "u_paral_stddev_resolution": r"$std dev \sigma({u_{\parallel}})$ [GeV]",
    "u_paral_stddev_scaled_resolution": r"$std dev \sigma({u_{\parallel}}) / (-<u_{\parallel}/q_{T}>)$",
}


qT_bins = np.array(
    [
        0,
        5,
        10,
        15,
        20,
        25,
        30,
        35,
        40,
        45,
        50,
        55,
        60,
        65,
        70,
        80,
        90,
        100,
        125,
        150,
        200,
        300,
        400,
        500,
    ]
)


hep.style.use("CMS")
color_dict = list(hep.style.CMS["axes.prop_cycle"])
color_list = [cycle["color"] for cycle in color_dict] + ["black", "darkgreen", "blue", "lightgreen"]
