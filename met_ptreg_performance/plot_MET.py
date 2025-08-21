import os
import sys
import re
from matplotlib import pyplot as plt
import matplotlib
import numpy as np
import awkward as ak
import argparse

import mplhep as hep

from hist import intervals
import hist
from hist import Hist



from utils.plot.get_columns_from_files import get_columns_from_files
from plot_config import total_var_dict,color_list  #  var_dict, color_list, ranges, log_dict


parser = argparse.ArgumentParser(description="Plot MET distributions from coffea files")
parser.add_argument(
    "-i",
    "--input-dir",
    type=str,
    required=True,
    help="Input directory for data with coffea files",
)
parser.add_argument("-o", "--output", type=str, help="Output directory", default="")

args = parser.parse_args()

outputdir = args.output if args.output else "plots_MET"
# Create output directory if it does not exist
if not os.path.exists(outputdir):
    os.makedirs(outputdir)


def plot_from_columns(cat_col):
    for category, col_var in cat_col.items():

        for quantity_name, var_dict in total_var_dict.items():

            fig, (ax, ax_ratio) = plt.subplots(
                2,
                1,
                figsize=[13, 13],
                sharex=True,
                gridspec_kw={"height_ratios": [2.5, 1]},
            )

            plot_name = var_dict["plot_name"]
            variables = var_dict["variables"]
            range_plot = var_dict["range"]
            log = var_dict["log"]

            for i, variable in enumerate(variables):
                print(f"Plotting {category} - {quantity_name} - {variable}")
                col_num = col_var[variable]
                col_den = col_var[variables[0]]
                var_name = (
                    variable.split("_")[0]
                    if "MuonGood" not in variable
                    else f'u {variable.split("_")[0]}'
                )

                mask_range = (col_den > range_plot[0]) & (col_den < range_plot[1])
                
                col_num = col_num[mask_range]
                col_den = col_den[mask_range]
                weight= col_var["weight"][mask_range]
                
                hist_den = Hist.new.Reg(30, range_plot[0], range_plot[1], name=var_name, flow=False).Weight()
                hist_num = Hist.new.Reg(30, range_plot[0], range_plot[1], name=var_name, flow=False).Weight()
                hist_den.fill(col_den, weight=weight)
                hist_num.fill(col_num, weight=weight)
                
                bins=hist_den.axes[0].edges
                bins_center = (hist_den.axes[0].edges[1:] + hist_den.axes[0].edges[:-1]) / 2
                
                ratio, err_ratio_up, err_ratio_down = hep.get_comparison(
                    hist_num, hist_den, comparison="split_ratio" if i == 0 else "ratio"
                )
                hep.histplot(
                    hist_num,
                    w2method="sqrt",
                    histtype="step",
                    label=var_name,
                    ax=ax,
                    color=color_list[i],
                )
                ax.set_xlabel("")

                if i == 0:
                    ax_ratio.axhline(y=1, linestyle="--",color=color_list[i],)
                    ax_ratio.fill_between(
                        bins_center,
                        1 - err_ratio_up, # ensure that the errors are symmetric
                        1 + err_ratio_up, # ensure that the errors are symmetric
                        # color="grey",
                        alpha=0.5,
                        color=color_list[i],
                    )
                else:
                    hep.histplot(
                        ratio,
                        bins=bins,
                        yerr=err_ratio_up,
                        histtype="errorbar",
                        label=var_name,
                        ax=ax_ratio,
                        color=color_list[i],
                    )
                    
                del col_den, col_num

            ax.legend(loc="upper right")
            ax.set_yscale("log" if log else "linear")
            hep.cms.lumitext(r"(13.6 TeV)", ax=ax)
            hep.cms.text(text="Simulation Preliminary", ax=ax)

            ax_ratio.set_xlabel(plot_name)
            ax.set_ylabel("Events")
            ax_ratio.set_ylabel(var_dict["ratio_label"])

            ax.grid()
            ax_ratio.grid()
            # ax_ratio.set_ylim(0.5, 1.5)
            ax_ratio.set_yscale("log" if log else "linear")
            ax.set_ylim(
                top=(1.7 * ax.get_ylim()[1] if not log else ax.get_ylim()[1] ** (1.7))
            )

            fig.savefig(
                os.path.join(outputdir, f"{category}_{quantity_name}.png"),
                bbox_inches="tight",
                dpi=300,
            )
            fig.savefig(
                os.path.join(outputdir, f"{category}_{quantity_name}.pdf"),
                bbox_inches="tight",
                dpi=300,
            )
            fig.savefig(
                os.path.join(outputdir, f"{category}_{quantity_name}.svg"),
                bbox_inches="tight",
                dpi=300,
            )
            plt.close(fig)


if __name__ == "__main__":

    inputfiles_data = [
        os.path.join(args.input_dir, file)
        for file in os.listdir(args.input_dir)
        if file.endswith(".coffea")
    ]

    cat_col, total_datasets_list = get_columns_from_files(inputfiles_data)
    print(f"Total datasets found: {total_datasets_list}")
    plot_from_columns(cat_col)
