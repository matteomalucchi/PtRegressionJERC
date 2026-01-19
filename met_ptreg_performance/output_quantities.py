import numpy as np
from pocket_coffea.parameters.histograms import (
    HistConf,
    Axis,
)
from pocket_coffea.lib.columns_manager import ColOut

bins_nPV = list(np.arange(0, 80, 1))
bins_qT = [
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
# save finely between -5 and 5, more coarsely outside
bins_response = list(
    np.concatenate(
        (np.linspace(-20, -5, 40), np.linspace(-4.999, 4.999, 10000), np.linspace(5, 20, 40))
    )
)
bins_u = list(
    np.concatenate(
        (
            np.linspace(-500, -200, 40),
            np.linspace(-199.999, 199.999, 10000),
            np.linspace(200, 500, 40),
        )
    )
)


def get_met_columns():
    met_vars = ["pt", "phi"]
    recoil_vars = ["pt", "phi", "u_perp_predict", "u_paral_predict", "response"]
    met_cols = []
    for recoil, vars_col in zip(["u", ""], [recoil_vars, met_vars]):
        # for raw in ["Raw", ""]:
        for raw in ["Raw"]:
            for type1 in [
                "",
                "-Type1",
                "-Type1JEC",
                "-Type1CorrMET",
                "-Type1CorrMETUncorrected",
                "-Type1PNetCorrMET",
                "-Type1PNetPlusNeutrinoCorrMET",
            ]:

                met_cols.append(ColOut(f"{recoil}{raw}PuppiMET{type1}", vars_col))

        met_cols.append(ColOut(f"{recoil}PuppiMET", vars_col))

    print("Total columns to be stored: ", met_cols)

    return met_cols


def get_met_variables():
    met_vars = {}
    for binning_variable in ["nPV", "qT"]:
        for variable in ["u_perp_predict", "u_paral_predict", "response"]:
            # define 2D histograms with binning_variable and the variable
            met_vars[f"{variable}_vs_{binning_variable}"] = HistConf(
                axes=[
                    Axis(
                        coll="PV" if binning_variable == "nPV" else "ll",
                        field="npvs" if binning_variable == "nPV" else "pt",
                        bins=bins_nPV if binning_variable == "nPV" else bins_qT,
                        label=(
                            "# PVs" if binning_variable == "nPV" else "Z $q_{T}$ [GeV]"
                        ),
                        pos=None,
                        type="variable",
                    ),
                    Axis(
                        coll="uRawPuppiMET",
                        field=variable,
                        bins=bins_response if variable == "response" else bins_u,
                        label=variable,
                        pos=None,
                        type="variable",
                    ),
                ],
            )

    # print("Total variables to be stored: ", met_vars)
    return met_vars

    # # define 2D histograms with PV_npvs, ll_pt and the mean accumulator of the variable
    # tot_vars[variable] = HistConf(
    #     axes=[
    #         Axis(
    #             coll="PV",
    #             field="npvs",
    #             bins=bins_nPV,
    #             label="# PVs",
    #             pos=None,
    #         ),
    #         Axis(
    #             coll="ll",
    #             field="pt",
    #             bins=bins_qT,
    #             label="Z $q_{T}$ [GeV]",
    #             pos=None,
    #         ),
    #     ],
    #     no_weights=False,
    #     accumulator_type="WeightedMean",
    # )
