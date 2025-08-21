import os
import sys
import re
from matplotlib import pyplot as plt
from coffea.util import load
from omegaconf import OmegaConf
import numpy as np
import awkward as ak
import mplhep as hep
import argparse
from collections import defaultdict
import functools
from hist import Hist
from multiprocessing import Pool
import matplotlib

from utils.plot.get_columns_from_files import get_columns_from_files
from plot_config import response_var_name_dict, qT_bins, color_list, met_list
from utils.plot.plotting import plot_1d_histograms

# hep.style.use("CMS")
# color_dict = list(hep.style.CMS["axes.prop_cycle"])
# color_list = [cycle["color"] for cycle in color_dict]

parser = argparse.ArgumentParser(description="Plot MET distributions from coffea files")
parser.add_argument(
    "-i",
    "--input-dir",
    type=str,
    required=True,
    help="Input directory for data with coffea files",
)
parser.add_argument(
    "--histo",
    action="store_true",
    help="If set, will plot 1d and 2d histograms of the recoil variables",
)
parser.add_argument(
    "-w",
    "--workers",
    type=int,
    default=1,
    help="Number of workers for multiprocessing (default: 1, no multiprocessing)",
)
parser.add_argument("-o", "--output", type=str, help="Output directory", default="")

args = parser.parse_args()


qT_bin_centers = []
for i in range(1, len(qT_bins)):
    qT_bin_centers.append((qT_bins[i] + qT_bins[i - 1]) / 2.0)
qT_bin_centers = np.array(qT_bin_centers)


outputdir = args.output if args.output else "plots_MET_response"
# Create output directory if it does not exist
if not os.path.exists(outputdir):
    os.makedirs(outputdir)

histograms_dir = os.path.join(outputdir, "1d_histograms")
if not os.path.exists(histograms_dir):
    os.makedirs(histograms_dir)
histograms_2d_dir = os.path.join(outputdir, "2d_histograms")
if not os.path.exists(histograms_2d_dir):
    os.makedirs(histograms_2d_dir)


def weighted_mean(x, w):
    mean_w = np.average(x, weights=w)
    # Effective number of entries
    n_eff = (np.sum(w)) ** 2 / np.sum(w**2)

    # Weighted variance
    variance_w = np.sum(w * (x - mean_w) ** 2) / np.sum(w)

    # Standard error
    sem_w = np.sqrt(variance_w / n_eff)

    return mean_w, sem_w


def weighted_std_dev(x, w):
    mean_w = np.average(x, weights=w)
    var_w = np.sum(w * (x - mean_w) ** 2) / np.sum(w)
    std_w = np.sqrt(var_w)

    n_eff = (np.sum(w) ** 2) / np.sum(w**2)
    std_err_w = std_w / np.sqrt(2 * (n_eff - 1))

    return std_w, std_err_w


def compute_u_info(u_i, weights_i, distribution_name, hists_dict):

    # compute mean
    mean_u_i, err_mean_u_i = weighted_mean(u_i, weights_i)

    hists_dict[f"{distribution_name}_mean"][0].append(mean_u_i)
    hists_dict[f"{distribution_name}_mean"][1].append(err_mean_u_i)

    # compute quantiles
    hists_dict[f"{distribution_name}_quantile_resolution"][0].append(
        (np.quantile(u_i, 0.84) - np.quantile(u_i, 0.16)) / 2.0,
    )
    # TODO: compute error on quantile
    hists_dict[f"{distribution_name}_quantile_resolution"][1].append(0)

    # compute standard deviation
    stddev_u_i, err_stddev_u_i = weighted_std_dev(u_i, weights_i)
    hists_dict[f"{distribution_name}_stddev_resolution"][0].append(stddev_u_i)
    hists_dict[f"{distribution_name}_stddev_resolution"][1].append(err_stddev_u_i)


def create_hist(hists_dict, qT_arr, u, weights, distribution_name, bin_edges):

    h = (
        Hist.new.Var(qT_bins, name="qT", label=r"Z q$_{\mathrm{T}}$ [GeV]", flow=False)
        .Var(bin_edges, name=distribution_name, label=distribution_name, flow=False)
        .Weight()
    )

    h.fill(
        qT_arr,
        u,
        weight=weights,
    )
    hists_dict[f"{distribution_name}"] = h


def create_reponses_info(qT_arr, u_dict, weights):
    # compute mean of all metrics in summary

    # max_x=200 # max qT value
    # x_n=20 #number of bins
    # bin_edges=np.arange(0, max_x, 10)

    bin_edges = qT_bins

    inds = np.digitize(qT_arr, bin_edges)

    all_responses = {}
    all_hists = {}
    for met_type in met_list:
        print(met_type)
        assert met_type in u_dict, f"MET type {met_type} not found in u_dict"
        
        for var_name in u_dict[met_type]:
            if "u_perp" in var_name:
                u_perp_arr = u_dict[met_type][var_name]
            elif "u_paral" in var_name:
                u_par_arr = u_dict[met_type][var_name]
            elif "response" in var_name:
                R_arr = u_dict[met_type][var_name]
        # u_perp_arr = u_dict[met_type][0]
        # u_par_arr = u_dict[met_type][1]
        # R_arr = u_dict[met_type][2]

        # print(R_arr, len(R_arr))
        # print(u_perp_arr, len(u_perp_arr))
        # print(u_par_arr, len(u_par_arr))

        # R_hist = []
        # u_perp_quantile_hist = []
        # u_perp_scaled_quantile_hist = []
        # u_perp_stddev_hist = []
        # u_perp_stddev_scaled_hist = []
        # u_par_quantile_hist = []
        # u_par_scaled_quantile_hist = []
        # u_par_stddev_hist = []
        # u_par_stddev_scaled_hist = []

        # err_R_hist = []
        # err_u_perp_quantile_hist = []
        # err_u_perp_scaled_quantile_hist = []
        # err_u_perp_stddev_hist = []
        # err_u_perp_stddev_scaled_hist = []
        # err_u_par_quantile_hist = []
        # err_u_par_scaled_quantile_hist = []
        # err_u_par_stddev_hist = []
        # err_u_par_stddev_scaled_hist = []

        # R_hist = defaultdict(list)
        # u_perp_hist = defaultdict(list)
        # u_perp_scaled_hist = defaultdict(list)
        # u_par_hist = defaultdict(list)
        # u_par_scaled_hist = defaultdict(list)

        all_hists[met_type] = {}
        R_bin_edges = np.linspace(-2, 2, 30)
        u_bin_edges = np.linspace(-200, 200, 30)

        create_hist(all_hists[met_type], qT_arr, R_arr, weights, "R", R_bin_edges)
        create_hist(
            all_hists[met_type], qT_arr, u_perp_arr, weights, "u_perp", u_bin_edges
        )
        # TODO: fix the scaling
        create_hist(
            all_hists[met_type],
            qT_arr,
            u_perp_arr / R_arr,
            weights,
            "u_perp_scaled",
            u_bin_edges,
        )
        create_hist(
            all_hists[met_type], qT_arr, u_par_arr, weights, "u_paral", u_bin_edges
        )
        # TODO: fix the scaling
        create_hist(
            all_hists[met_type],
            qT_arr,
            u_par_arr / R_arr,
            weights,
            "u_paral_scaled",
            u_bin_edges,
        )
        # breakpoint()

        all_responses[met_type] = defaultdict(lambda: [[], []])
        # continue

        for i in range(1, len(bin_edges)):
            weights_i = weights[np.where(inds == i)[0]]

            # Response
            R_i = R_arr[np.where(inds == i)[0]]
            av_R_i, _ = weighted_mean(R_i, weights_i)
            compute_u_info(R_i, weights_i, "R", all_responses[met_type])
            # compute mean and standard deviation
            # R_hist.append(av_R_i)
            # err_R_hist.append(err_R_i)

            # U perpendicular
            u_perp_i = u_perp_arr[np.where(inds == i)[0]]
            compute_u_info(
                u_perp_i,
                weights_i,
                "u_perp",
                all_responses[met_type],
            )

            u_perp_scaled_i = u_perp_i / av_R_i
            compute_u_info(
                u_perp_scaled_i,
                weights_i,
                "u_perp_scaled",
                all_responses[met_type],
            )

            # U parallel
            u_par_i = u_par_arr[np.where(inds == i)[0]]
            compute_u_info(
                u_par_i,
                weights_i,
                "u_par",
                all_responses[met_type],
            )

            u_par_scaled_i = u_par_i / av_R_i
            compute_u_info(
                u_par_scaled_i,
                weights_i,
                "u_par_scaled",
                all_responses[met_type],
            )

            # compute_u_hists(
            #     u_perp_i,
            #     weights_i,
            #     u_perp_quantile_hist,
            #     err_u_perp_quantile_hist,
            #     u_perp_stddev_hist,
            #     err_u_perp_stddev_hist,
            # )
            # compute_u_hists(
            #     u_perp_scaled_i,
            #     weights_i,
            #     u_perp_scaled_quantile_hist,
            #     err_u_perp_scaled_quantile_hist,
            #     u_perp_stddev_scaled_hist,
            #     err_u_perp_stddev_scaled_hist,
            # )

            # ## compute quantiles
            # u_perp_quantile_hist.append(
            #     (np.quantile(u_perp_i, 0.84) - np.quantile(u_perp_i, 0.16)) / 2.0
            # )
            # err_u_perp_quantile_hist.append(0)
            # u_perp_scaled_quantile_hist.append(
            #     (
            #         np.quantile(u_perp_scaled_i, 0.84)
            #         - np.quantile(u_perp_scaled_i, 0.16)
            #     )
            #     / 2.0
            # )
            # err_u_perp_scaled_quantile_hist.append(0)

            # ## compute standard deviation
            # stdev_u_perp_i, err_u_perp_stddev_i = weighted_std_dev(u_perp_i, weights_i)
            # u_perp_stddev_hist.append(stdev_u_perp_i)
            # err_u_perp_stddev_hist.append(err_u_perp_stddev_i)
            # stdev_u_perp_scaled_i, err_u_perp_stddev_scaled_i = weighted_std_dev(
            #     u_perp_scaled_i, weights_i
            # )
            # u_perp_stddev_scaled_hist.append(stdev_u_perp_scaled_i)
            # err_u_perp_stddev_scaled_hist.append(err_u_perp_stddev_scaled_i)

            # compute_u_hists(
            #     u_par_i,
            #     weights_i,
            #     u_par_quantile_hist,
            #     err_u_par_quantile_hist,
            #     u_par_stddev_hist,
            #     err_u_par_stddev_hist,
            # )
            # compute_u_hists(
            #     u_par_scaled_i,
            #     weights_i,
            #     u_par_scaled_quantile_hist,
            #     err_u_par_scaled_quantile_hist,
            #     u_par_stddev_scaled_hist,
            #     err_u_par_stddev_scaled_hist,
            # )

            # u_par_quantile_hist.append(
            #     (np.quantile(u_par_i, 0.84) - np.quantile(u_par_i, 0.16)) / 2.0
            # )
            # u_par_scaled_quantile_hist.append(
            #     (np.quantile(u_par_scaled_i, 0.84) - np.quantile(u_par_scaled_i, 0.16))
            #     / 2.0
            # )
            # err_u_par_quantile_hist.append(0)
            # err_u_par_scaled_quantile_hist.append(0)

        # u_perp_quantile_resolution = np.histogram(qT_bin_centers, bins=bin_edges, weights=u_perp_quantile_hist)
        # u_perp_scaled_quantile_resolution = np.histogram(
        #     qT_bin_centers, bins=bin_edges, weights=u_perp_scaled_quantile_hist
        # )
        # u_par_quantile_resolution = np.histogram(qT_bin_centers, bins=bin_edges, weights=u_par_quantile_hist)
        # u_par_scaled_quantile_resolution = np.histogram(
        #     qT_bin_centers, bins=bin_edges, weights=u_par_scaled_quantile_hist
        # )
        # R = np.histogram(qT_bin_centers, bins=bin_edges, weights=R_hist)

        # all_hists[key] = {
        #     "u_perp_quantile_resolution": u_perp_quantile_resolution,
        #     "u_perp_scaled_quantile_resolution": u_perp_scaled_quantile_resolution,
        #     "u_par_quantile_resolution": u_par_quantile_resolution,
        #     "u_par_scaled_quantile_resolution": u_par_scaled_quantile_resolution,
        #     "R": R,
        # }

        # all_hists[key] = {
        #     "u_perp_quantile_resolution": (
        #         u_perp_quantile_hist,
        #         err_u_perp_quantile_hist,
        #     ),
        #     "u_perp_scaled_quantile_resolution": (
        #         u_perp_scaled_quantile_hist,
        #         err_u_perp_scaled_quantile_hist,
        #     ),
        #     "u_perp_stddev_resolution": (u_perp_stddev_hist, err_u_perp_stddev_hist),
        #     "u_perp_stddev_scaled_resolution": (
        #         u_perp_stddev_scaled_hist,
        #         err_u_perp_stddev_scaled_hist,
        #     ),
        #     "u_par_quantile_resolution": (u_par_quantile_hist, err_u_par_quantile_hist),
        #     "u_par_scaled_quantile_resolution": (
        #         u_par_scaled_quantile_hist,
        #         err_u_par_scaled_quantile_hist,
        #     ),
        #     "u_par_stddev_resolution": (u_par_stddev_hist, err_u_par_stddev_hist),
        #     "u_par_stddev_scaled_resolution": (
        #         u_par_stddev_scaled_hist,
        #         err_u_par_stddev_scaled_hist,
        #     ),
        #     "R": (R_hist, err_R_hist),
        # }

    # change the gerarchy of the keys
    reponses_dict = {}
    for met_type in all_responses:
        for var_name in all_responses[met_type]:
            if var_name not in reponses_dict:
                reponses_dict[var_name] = {}
            reponses_dict[var_name][met_type] = all_responses[met_type][var_name]

    hists_dict = {}
    for met_type in all_hists:
        for var_name in all_hists[met_type]:
            if var_name not in hists_dict:
                hists_dict[var_name] = {}
            hists_dict[var_name][met_type] = all_hists[met_type][var_name]

    # breakpoint()

    return reponses_dict, hists_dict


def plot_reponses(reponses_dict, cat):
    for var_name in reponses_dict:
        print(f"Plotting response for {var_name} in category {cat}")
        # [
        #     "u_perp_quantile_resolution",
        #     "u_perp_scaled_quantile_resolution",
        #     "u_par_quantile_resolution",
        #     "u_par_scaled_quantile_resolution",
        #     "R",
        # ]:
        fig, ax = plt.subplots()
        for i, met_type in enumerate(reponses_dict[var_name]):
            ax.errorbar(
                qT_bin_centers,
                reponses_dict[var_name][met_type][0],
                xerr=(qT_bins[1:] - qT_bins[:-1]) / 2.0,
                yerr=reponses_dict[var_name][met_type][1],
                label=met_type,
                color=color_list[i],
                fmt=".",
            )
        ax.legend(loc="best")
        ax.set_xlabel(r"Z q$_{\mathrm{T}}$ [GeV]")
        ax.set_ylabel(
            var_name
            if var_name not in response_var_name_dict
            else response_var_name_dict[var_name]
        )
        hep.cms.lumitext(r"(13.6 TeV)", ax=ax)
        hep.cms.text(text="Preliminary", ax=ax)
        fig.savefig(f"{outputdir}/{cat}_{var_name}.png", bbox_inches="tight", dpi=300)
        fig.savefig(f"{outputdir}/{cat}_{var_name}.pdf", bbox_inches="tight", dpi=300)
        fig.savefig(f"{outputdir}/{cat}_{var_name}.svg", bbox_inches="tight", dpi=300)
        plt.close(fig)


# def plot_2d_histogram


def plot_2d_response_histograms(hists_dict, cat):
    for var_name in hists_dict:
        print(f"Plotting 2d histogram for {var_name} in category {cat}")
        for met_type in hists_dict[var_name]:
            fig, ax = plt.subplots()
            hist = hists_dict[var_name][met_type]

            hep.hist2dplot(
                hist,
                ax=ax,
                # xaxis="qT",
                # yaxis=var_name,
                label=met_type,
                cmap="viridis",
                norm=matplotlib.colors.LogNorm(),
            )
            ax.set_xlabel(r"Z q$_{\mathrm{T}}$ [GeV]")
            ax.set_ylabel(
                var_name
                if var_name not in response_var_name_dict
                else response_var_name_dict[var_name]
            )
            ax.legend(loc="best")
            hep.cms.lumitext(r"(13.6 TeV)", ax=ax)
            hep.cms.text(text="Preliminary", ax=ax)
            fig.savefig(
                f"{histograms_2d_dir}/2d_histo_{cat}_{var_name}_{met_type}.png",
                bbox_inches="tight",
                dpi=300,
            )
            fig.savefig(
                f"{histograms_2d_dir}/2d_histo_{cat}_{var_name}_{met_type}.pdf",
                bbox_inches="tight",
                dpi=300,
            )
            fig.savefig(
                f"{histograms_2d_dir}/2d_histo_{cat}_{var_name}_{met_type}.svg",
                bbox_inches="tight",
                dpi=300,
            )
            plt.close(fig)


# def old_plot_1d_histograms(hists_dict, cat):
#     # for each bin on qT, plot the distribution of the variable
#     for i in range(len(qT_bins) - 1):
#         bin_edges_string = f"{qT_bins[i]}_{qT_bins[i+1]}"
#         for var_name in hists_dict:
#             fig, ax = plt.subplots()
#             for met_type in hists_dict[var_name]:
#                 print(
#                     f"Plotting 1d histogram for {var_name} in {met_type} for qT bin {bin_edges_string}"
#                 )
#                 hist = hists_dict[var_name][met_type]
#                 # Select the bin corresponding to the current qT bin
#                 hist_1d_u = hist[{"qT": i}]
#                 hep.histplot(
#                     hist_1d_u,
#                     ax=ax,
#                     label=met_type,
#                     # color=hep.style.CMS.colors[key],
#                     histtype="step",
#                 )
#             ax.legend(loc="best")
#             ax.set_xlabel(
#                 var_name
#                 if var_name not in response_var_name_dict
#                 else response_var_name_dict[var_name]
#             )
#             ax.set_ylabel("Events")
#             ax.set_ylim(top=1.7 * ax.get_ylim()[1])
#             hep.cms.lumitext(r"(13.6 TeV)", ax=ax)
#             hep.cms.text(text="Preliminary", ax=ax)
#             fig.savefig(
#                 f"{outputdir}/{cat}_{var_name}_{bin_edges_string}.png",
#                 bbox_inches="tight",
#                 dpi=300,
#             )
#             fig.savefig(
#                 f"{outputdir}/{cat}_{var_name}_{bin_edges_string}.pdf",
#                 bbox_inches="tight",
#                 dpi=300,
#             )
#             fig.savefig(
#                 f"{outputdir}/{cat}_{var_name}_{bin_edges_string}.svg",
#                 bbox_inches="tight",
#                 dpi=300,
#             )
#             plt.close(fig)


def plot_1d_histograms_parallel(plotting_info, log_scale, ratio_label):
    hists_dict, output_name, var_label = plotting_info
    print(f"Plotting 1d histogram in parallel for {output_name} with label {var_label}")
    plot_1d_histograms(
        hists_dict=hists_dict,
        output_name=output_name,
        var_label=var_label,
        log_scale=log_scale,
        ratio_label=ratio_label,
    )


def plot_1d_response_histograms(hists_dict, cat):
    # for each bin on qT, plot the distribution of the variable
    plotting_info_list = []
    for i in range(len(qT_bins) - 1):
        bin_edges_string = f"{qT_bins[i]}_{qT_bins[i+1]}"
        for var_name in hists_dict:
            hist_1d_dict = {}
            for j, met_type in enumerate(hists_dict[var_name]):
                hist = hists_dict[var_name][met_type]
                # Select the bin corresponding to the current qT bin
                hist_1d_u = hist[{"qT": i}]
                ratio_hist_den = True if met_type == "RawPuppiMET" else False
                hist_1d_dict[met_type] = (hist_1d_u, ratio_hist_den, color_list[j])
                output_name = f"{histograms_dir}/{cat}_{var_name}_{bin_edges_string}"
                var_label = (
                    var_name
                    if var_name not in response_var_name_dict
                    else response_var_name_dict[var_name]
                )
                if args.workers > 1:
                    plotting_info_list.append((hist_1d_dict, output_name, var_label))
                else:
                    plot_1d_histograms(
                        hists_dict=hist_1d_dict,
                        output_name=output_name,
                        var_label=var_label,
                        log_scale=False,
                        ratio_label=None,
                    )

    if args.workers > 1:
        with Pool(args.workers) as pool:
            pool.starmap(
                functools.partial(
                    plot_1d_histograms, log_scale=False, ratio_label=None
                ),
                plotting_info_list,
            )


def make_plots(cat_col):
    for cat in cat_col:
        print(f"Processing category: {cat}")
        col_dict = cat_col[cat]
        v_qT = col_dict["ll_pt"]

        u_dict = {}
        for var in col_dict:
            if "_MuonGood" in var and any(
                x in var for x in ["u_perp_predict", "u_paral_predict", "response"]
            ):
                coll = var.split("_")[0]
                if coll not in u_dict:
                    u_dict[coll] = {}
                print(var)
                u_dict[coll][var] = col_dict[var]
            elif "weight" in var:
                weights = col_dict[var]
        # breakpoint()
        reponses_dict, hists_dict = create_reponses_info(v_qT, u_dict, weights)
        plot_reponses(reponses_dict, cat)
        if args.histo:
            plot_2d_response_histograms(hists_dict, cat)
            plot_1d_response_histograms(hists_dict, cat)


if __name__ == "__main__":

    inputfiles_data = [
        os.path.join(args.input_dir, file)
        for file in os.listdir(args.input_dir)
        if file.endswith(".coffea")
    ]

    cat_col, total_datasets_list = get_columns_from_files(inputfiles_data)
    print(f"Total datasets found: {total_datasets_list}")
    print(cat_col)

    make_plots(cat_col)
