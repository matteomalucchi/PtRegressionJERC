from coffea.util import load
from utils_configs.plot.HEPPlotter import HEPPlotter
import copy
import os
import argparse


FILE_DICT={
    "2022_preEE":"/work/mmalucch/out_MET/out_option_6_DYto2L-4Jets_MLL-50-v15_2022_preEE/plot_inclusive_distributions/histograms_baseline.coffea",
    "2022_postEE":"/work/mmalucch/out_MET/out_option_6_DYto2L-4Jets_MLL-50-v15_2022_postEE/plot_inclusive_distributions/histograms_baseline.coffea",
    "2023_preBPix":"/work/mmalucch/out_MET/out_option_6_DYto2L-4Jets_MLL-50-v15_2023_preBPix/plot_inclusive_distributions/histograms_baseline.coffea",
    "2023_postBPix":"/work/mmalucch/out_MET/out_option_6_DYto2L-4Jets_MLL-50-v15_2023_postBPix/plot_inclusive_distributions/histograms_baseline.coffea",
    "2024":"/work/mmalucch/out_MET/out_option_6_DYto2L-4Jets_MLL-50-v15_2024/plot_inclusive_distributions/histograms_baseline.coffea",
}


parser = argparse.ArgumentParser(description="Compare histograms across years.")
parser.add_argument(
    "-o",
    "--output-dir",
    type=str,
    required=True,
    help="Output directory for plots.",
)
parser.add_argument(
    "-v",
    "--var-names",
    nargs="+",
    default=["ll_pt_RegBin", "PV_npvs_RegBin"],
    help="Variable names to plot.",
)
args=parser.parse_args()


def main():
    # Load all histograms
    loaded = {}
    for year, file_name in FILE_DICT.items():
        d = load(file_name)
        loaded[year] = d["inclusive_histos"]

    years = list(FILE_DICT.keys())
    reference_year = years[0]

    plot_dict = {}
    for var_name in args.var_names:
        # Use the first year as the base structure
        plot_dict[var_name] = copy.deepcopy(loaded[reference_year][var_name])
        plot_dict[var_name]["series_dict"].pop("bin_var")
        plot_dict[var_name]["legend"] = True

        for i, year in enumerate(years):
            plot_dict[var_name]["series_dict"][year] = copy.deepcopy(
                loaded[year][var_name]["series_dict"][list(loaded[year][var_name]["series_dict"].keys())[0]]
            )
            plot_dict[var_name]["series_dict"][year]["style"]["legend_name"] = year
            plot_dict[var_name]["series_dict"][year]["style"]["is_reference"] = (i == 0)

    os.makedirs(args.output_dir, exist_ok=True)

    for var_name, info in plot_dict.items():
        print(f"Plotting {var_name}")
        p = (
            HEPPlotter()
            .set_plot_config(figsize=(12, 12), lumitext="(13.6 TeV)")
            .set_output(os.path.join(args.output_dir, var_name))
            .set_labels(info["xlabel"], "a.u.")
            .set_extra_kwargs(density=True)
            .set_options(
                y_log=info["y_log"],
                set_ylim=True,
                ylim_top_value=10,
                ylim_bottom_value=1e-7,
                y_log_ratio=True if "PV" in var_name else False,
                ylim_ratio_top_value=100 if "PV" in var_name else None,
                ylim_ratio_bottom_value=0.01 if "PV" in var_name else None,
                legend=info.get("legend", True),
                legend_loc="upper left",
                legend_font_size=16,
                normalize_1d_histo=True,
            )
            .set_data(info["series_dict"], plot_type="1d")
            .add_mean_std(x=0.95, y=0.95, fontsize=16)
        )
        p.run()


if __name__ == "__main__":
    main()