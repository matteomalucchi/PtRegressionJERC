import logging

logger = logging.getLogger("matplotlib")
logger.setLevel(logging.WARNING)  # suppress INFO
logger.propagate = False

import os
import numpy as np
from multiprocessing import Pool
import argparse
import hist
from coffea.util import load
from scipy.optimize import curve_fit
from scipy.stats import chi2 as chi2_dist

from utils_configs.plot.HEPPlotter import HEPPlotter
import met_ptreg_performance.helpers as helpers

from mc_truth_ptreg_jerc.response_plot.plot_config_jer_mc import (
    BIN_VARIABLES,
    BIN_VARIABLES_NEUTRINO,
    MAPPING_VARIABLES,
    RESPONSE_VARIABLES,
    RESPONSE_VARIABLES_NEUTRINO,
)

from mc_truth_ptreg_jerc.response_plot.confidence import Confidence_numpy


parser = argparse.ArgumentParser(
    description="Plot jet resolution from pre-computed histograms"
)
parser.add_argument(
    "-i",
    "--input-file",
    type=str,
    required=True,
    help="Input coffea file with pre-computed histograms",
)
parser.add_argument(
    "-w",
    "--workers",
    type=int,
    default=1,
    help="Number of workers for multiprocessing (default: 1, no multiprocessing)",
)
parser.add_argument("-o", "--output", type=str, help="Output directory", default="./")

args = parser.parse_args()


def plot_resolution_vs_x_variable(
    response_types_dict, response_vars, bin_var_configs, output_dir="", year=""
):
    """
    Plot resolution as a function of the variable marked with resolution_x_variable=True,
    grouped by bins of all other variables (resolution_x_variable=False).
    Overlays multiple response types on the same plot using HEPPlotter.

    Parameters
    ----------
    response_types_dict : dict
        Dictionary mapping response variable names to their resolution_result dicts
        (output from compute_binned_resolution_from_hists).
    response_vars : dict
        Dictionary with configuration for response variables.
    bin_var_configs : dict
        Configuration dict for binning variables.
    output_dir : str
        Output directory for saving plots (optional).
    year : str
        Year string for plot annotation.

    Returns
    -------
    dict
        Dictionary mapping each bin combination of non-x-axis variables to HEPPlotter objects.
    """

    if not response_types_dict:
        raise ValueError("response_types_dict cannot be empty")

    # Get the first resolution_result to extract bin structure
    first_response_var = list(response_types_dict.keys())[0]
    first_result = response_types_dict[first_response_var]

    bin_var_names = first_result["bin_var_names"]
    bin_edges_dict = first_result["bin_edges"]

    # Identify x-axis variable (resolution_x_variable=True)
    x_var = None
    y_vars = []  # Variables with resolution_x_variable=False
    x_var_idx = None

    for i, var_name in enumerate(bin_var_names):
        if bin_var_configs[var_name].get("resolution_x_variable", False):
            x_var = var_name
            x_var_idx = i
        else:
            y_vars.append((var_name, i))

    if x_var is None:
        raise ValueError("No variable with resolution_x_variable=True found")

    if len(y_vars) == 0:
        raise ValueError("No variables with resolution_x_variable=False found")

    # Get bin edges for x variable
    x_bin_edges = bin_edges_dict[x_var]
    x_bin_centers = (x_bin_edges[:-1] + x_bin_edges[1:]) / 2
    n_x_bins = len(x_bin_centers)

    # Determine shape of y-variable bins
    y_var_names = [v[0] for v in y_vars]
    y_bin_shape = tuple(len(bin_edges_dict[v]) - 1 for v in y_var_names)

    # Create plotters for each combination of y-variable bins
    plotters = {}

    for y_bin_idx in np.ndindex(y_bin_shape):
        # Collect data for all response types
        graph_data = {}

        for response_var, resolution_result in response_types_dict.items():
            resolutions = resolution_result["resolutions"]

            # Extract resolutions for this combination of y-variables across all x-bins
            resolutions_for_plot = []
            errors_for_plot = []
            valid_x_indices = []

            for x_idx in range(n_x_bins):
                # Build full bin index
                full_bin_idx = list(range(len(bin_var_names)))
                for i, (y_var_name, y_var_idx) in enumerate(y_vars):
                    full_bin_idx[y_var_idx] = y_bin_idx[i]
                full_bin_idx[x_var_idx] = x_idx

                full_bin_idx = tuple(full_bin_idx)

                # Get resolution for this bin
                if (
                    full_bin_idx in resolutions
                    and resolutions[full_bin_idx] is not None
                ):
                    resolutions_for_plot.append(resolutions[full_bin_idx])
                    errors_for_plot.append(0)
                    valid_x_indices.append(x_idx)

            if len(resolutions_for_plot) == 0:
                continue

            # Get x-values for valid bins
            x_values = x_bin_centers[valid_x_indices]

            # Add to graph data
            graph_data[response_var] = {
                "data": {
                    "x": [x_values, np.zeros_like(x_values)],
                    "y": [np.array(resolutions_for_plot), np.array(errors_for_plot)],
                },
                "style": {
                    "fmt": "o",
                    "color": response_vars[response_var].get("color"),
                    "markersize": 8,
                    "linewidth": 2,
                    "legend_name": response_vars[response_var].get(
                        "legend_name", response_var
                    ),
                },
            }

        if not graph_data:
            continue

        # Build filename with y-variable bin ranges
        filename_parts = ["resolution"]

        for y_var_name, y_idx in y_vars:
            bin_idx = y_bin_idx[y_vars.index((y_var_name, y_idx))]
            low_edge_str = f"{bin_edges_dict[y_var_name][bin_idx]}".replace(
                ".", "p"
            ).replace("-", "m")
            high_edge_str = f"{bin_edges_dict[y_var_name][bin_idx+1]}".replace(
                ".", "p"
            ).replace("-", "m")
            filename_parts.append(
                f"{bin_var_configs[y_var_name]['name_plot']}_{low_edge_str}to{high_edge_str}"
            )

        # Create HEPPlotter
        output_name = "_".join(filename_parts) if output_dir else None

        # Build annotation text from bin ranges
        annotation_text = r"anti-$k_{T}$ R=0.4 (PUPPI)"
        for y_var_name, y_idx in y_vars:
            bin_idx = y_bin_idx[y_vars.index((y_var_name, y_idx))]
            low_edge = bin_edges_dict[y_var_name][bin_idx]
            high_edge = bin_edges_dict[y_var_name][bin_idx + 1]
            label = bin_var_configs[y_var_name].get("label", y_var_name)
            annotation_text += f"\n{low_edge} < {label} < {high_edge}"

        plotter = (
            HEPPlotter(debug=True)
            .set_plot_config(cmstext_loc=2, lumitext=f"{year} (13.6 TeV)", cmstext_font_size=35)
            .set_labels(
                xlabel=f"{bin_var_configs[x_var]['label']}",
                ylabel="Jet Energy Resolution",
            )
            .set_data(graph_data, plot_type="graph")
            .set_options(grid=True, legend=True, legend_loc="top right")
            .set_output(f"{output_dir}/{output_name}")
            .add_annotation(
                0.05,
                0.78,
                annotation_text,
                horizontalalignment="left",
                verticalalignment="top",
            )
        )

        plotters[y_bin_idx] = plotter

    return plotters


def perform_linear_fit(h_mean_hist):
    """
    Perform a weighted linear fit of mean rho versus pileup.

    Parameters
    ----------
    h_mean_hist : hist.Hist
        A 1D histogram with Mean storage for rho versus Pileup_nPU.

    Returns
    -------
    dict
        Dictionary containing fit parameters, uncertainties, fit quality metrics,
        and a callable prediction function.
    """

    if h_mean_hist.ndim != 1:
        raise ValueError("perform_linear_fit requires a 1D Mean histogram")

    x_axis = h_mean_hist.axes[0]
    if not hasattr(x_axis, "centers"):
        raise ValueError("Histogram axis does not expose centers")

    x = np.array(x_axis.centers)
    view = h_mean_hist.view()
    y = np.array(view.value)
    variance = np.array(view.variance)

    if y.shape != x.shape:
        raise ValueError("Histogram x and y arrays must have the same shape")

    y_err = np.sqrt(variance)
    valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(y_err) & (y_err > 0)
    if not np.any(valid):
        raise ValueError("Histogram contains no valid points for fitting")

    x = x[valid]
    y = y[valid]
    y_err = y_err[valid]

    def linear_model(x, m, b):
        return m * x + b

    popt, pcov = curve_fit(linear_model, x, y, sigma=y_err, absolute_sigma=True)
    slope, intercept = popt
    slope_err, intercept_err = np.sqrt(np.diag(pcov))

    y_pred = linear_model(x, slope, intercept)
    residuals = y - y_pred

    weights = 1.0 / y_err**2
    ss_res = np.sum(weights * residuals**2)
    y_mean = np.sum(weights * y) / np.sum(weights)
    ss_tot = np.sum(weights * (y - y_mean) ** 2)

    r_squared = 1.0 - ss_res / ss_tot if ss_tot != 0 else np.nan
    chi2_stat = ss_res
    dof = max(len(x) - 2, 1)
    p_value = chi2_dist.sf(chi2_stat, dof)

    return {
        "slope": slope,
        "slope_err": slope_err,
        "intercept": intercept,
        "intercept_err": intercept_err,
        "r_squared": r_squared,
        "chi2": chi2_stat,
        "p_value": p_value,
        "n_points": len(x),
        "fit_func": lambda x_val: linear_model(x_val, slope, intercept),
    }


def _compute_resolution_from_hist(args):
    """
    Helper function for parallel processing of a single response variable histogram.
    Computes binned resolution from pre-computed hist objects.

    Parameters
    ----------
    args : tuple
        (response_var_name, hist_obj, bin_var_names, bin_var_configs, response_vars)

    Returns
    -------
    tuple
        (response_var_name, resolution_result_dict)
    """
    response_var_name, hist_obj, bin_var_names, bin_var_configs, response_vars = args

    n_response_bins = response_vars.get("N_bins", 50)

    if len(bin_var_names) == 0:
        raise ValueError("At least one bin variable is required")

    # Extract bin edges from config
    bin_edges_dict = {
        var_name: bin_var_configs[var_name]["bin_edges"]
        for var_name in bin_var_names
    }

    print(f"Processing resolution from histogram '{response_var_name}'...")

    # Extract bin shape
    bin_shape = tuple(len(edges) - 1 for edges in bin_edges_dict.values())

    # Initialize resolution containers
    resolutions = {}
    resolution_grid = np.full(bin_shape, np.nan)

    # Iterate over all bin combinations in the histogram
    for bin_idx in np.ndindex(bin_shape):
        try:
            # Build the slice specification for the histogram
            # hist objects support integer indexing for bins
            # bin_idx gives us the bin indices we want to slice
            sliced_hist = hist_obj[bin_idx]

            # The sliced histogram should be 1D (just the response axis)
            # Get the values from this 1D histogram
            counts = sliced_hist.view()

            if np.sum(counts) > 5:  # Only compute if sufficient data points
                # Extract axis info from sliced histogram
                response_axis = sliced_hist.axes[0]
                bin_centers = np.array(response_axis.centers)
                bin_width = response_axis.widths[0] if len(response_axis.widths) > 0 else 1.0

                resolution = Confidence_numpy(
                    hist=counts,
                    bins_mid=bin_centers,
                    bin_width=bin_width,
                    confLevel=0.68,
                )
                resolutions[bin_idx] = resolution
                resolution_grid[bin_idx] = resolution
        except Exception as e:
            print(f"Warning: Failed to extract resolution for bin {bin_idx}: {e}")

    return response_var_name, {
        "response_var": response_var_name,
        "bin_edges": bin_edges_dict,
        "bin_var_names": bin_var_names,
        "resolutions": resolutions,
        "resolution_grid": resolution_grid,
    }


def compute_binned_resolution_from_hists(
    hist_dict, bin_var_names, bin_var_configs, response_vars, n_workers=4
):
    """
    Compute resolution (width) from pre-computed histograms.

    Parameters
    ----------
    hist_dict : dict
        Dictionary mapping response_var_name -> hist.Hist object.
        Each histogram should have dimensions (bin_vars..., response_var).
    bin_var_names : list
        List of bin variable names.
    bin_var_configs : dict
        Configuration dict for binning variables.
    response_vars : dict
        Configuration dict for response variables.
    n_workers : int
        Number of parallel workers (default: 4).

    Returns
    -------
    dict
        Maps response_var_name -> resolution_result dict.
    """

    # Prepare arguments for parallel processing
    args_list = [
        (response_var_name, hist_obj, bin_var_names, bin_var_configs, response_vars[response_var_name])
        for response_var_name, hist_obj in hist_dict.items()
    ]

    # Use multiprocessing Pool for parallel computation
    response_tot_dict = {}

    if len(args_list) == 1 or n_workers == 1:
        # Single worker - no need for Pool overhead
        for args in args_list:
            var_name, result = _compute_resolution_from_hist(args)
            response_tot_dict[var_name] = result
    else:
        # Multi-worker processing
        with Pool(n_workers) as pool:
            results = pool.map(_compute_resolution_from_hist, args_list)
            for var_name, result in results:
                response_tot_dict[var_name] = result

    return response_tot_dict


def plot_rho_vs_pu(h_mean_rho, year, category):
    """
    Plot and fit rho vs pileup from pre-computed histogram.

    Parameters
    ----------
    h_mean_rho : hist.Hist
        1D histogram with Mean storage: (PU,) -> mean rho values.
    year : str
        Year string for annotation.
    category : str
        Data category string.

    Returns
    -------
    dict
        Fit results dictionary.
    """

    fit_results = perform_linear_fit(h_mean_rho)
    print(fit_results)

    x_axis = h_mean_rho.axes[0]
    mean_rho_vs_mu_dict = {
        "data": {
            "data": {
                "x": [(x_axis.edges[:-1] + x_axis.edges[1:]) / 2, x_axis.widths / 2],
                "y": [
                    h_mean_rho.view().value,
                    np.sqrt(h_mean_rho.view().variance),
                ],
            },
            "style": {
                "fmt": "o",
            },
        },
        "linear fit": {
            "data": {
                "x": [
                    np.linspace(x_axis.edges[0], x_axis.edges[-1], 100),
                    np.zeros(100),
                ],
                "y": [
                    fit_results["fit_func"](
                        np.linspace(x_axis.edges[0], x_axis.edges[-1], 100)
                    ),
                    np.zeros(100),
                ],
            },
            "style": {
                "linestyle": "-",
                "fmt": "",
            },
        },
    }

    p = (
        HEPPlotter()
        .set_plot_config(lumitext=f"{year} (13.6 TeV)")
        .set_options(
            set_ylim=False,
        )
        .set_output(f"{args.output}/rho_vs_pu_fit_{category}")
        .set_data(mean_rho_vs_mu_dict, plot_type="graph")
        .set_labels("$\\mu$", r"$< \rho >$")
        .add_annotation(
            x=0.05,
            y=0.7,
            s=f"$\\chi^2/ndf$ {fit_results['chi2']/(fit_results['n_points'] - 2):.3f}\n"
            + f"$p$-value {fit_results['p_value']:.3f}\n",
        )
    )
    p.run()

    return fit_results


def main():
    """
    Main function: load pre-computed histograms and produce resolution plots.
    """
    # Create output directory
    os.makedirs(args.output, exist_ok=True)

    print(f"Loading pre-computed histograms from {args.input_file}...")
    histograms = load(args.input_file)

    # Expected histogram structure in coffea file:
    # {
    #     'category_name': {
    #         'h_rho_pu': hist.Hist (PU,) with Mean storage,
    #         'h_ResponseJEC': hist.Hist (PU, eta, pt, ResponseJEC),
    #         'h_ResponseRaw': hist.Hist (PU, eta, pt, ResponseRaw),
    #         'h_ResponsePNetReg': hist.Hist (PU, eta, pt, ResponsePNetReg),
    #         'h_ResponseUparTReg': hist.Hist (PU, eta, pt, ResponseUparTReg),
    #         'h_ResponseJEC_neutrino': hist.Hist (PU, eta, pt, ResponseJEC_neutrino),
    #         'h_ResponseRaw_neutrino': hist.Hist (PU, eta, pt, ResponseRaw_neutrino),
    #         'h_ResponsePNetRegNeutrino': hist.Hist (PU, eta, pt, ResponsePNetRegNeutrino),
    #         'h_ResponseUparTRegNeutrino': hist.Hist (PU, eta, pt, ResponseUparTRegNeutrino),
    #     },
    # }

    # Mapping from RESPONSE_VARIABLES keys to histogram names
    response_to_hist_name = {
        "MatchedJets_ResponseJEC": "h_ResponseJEC",
        "MatchedJets_ResponseRaw": "h_ResponseRaw",
        "MatchedJets_ResponsePNetReg": "h_ResponsePNetReg",
        "MatchedJets_ResponseUparTReg": "h_ResponseUparTReg",
        "MatchedJetsNeutrino_ResponseJEC": "h_ResponseJEC_neutrino",
        "MatchedJetsNeutrino_ResponseRaw": "h_ResponseRaw_neutrino",
        "MatchedJetsNeutrino_ResponsePNetRegNeutrino": "h_ResponsePNetRegNeutrino",
        "MatchedJetsNeutrino_ResponseUparTRegNeutrino": "h_ResponseUparTRegNeutrino",
    }

    for category, hist_data in histograms.items():
        print(f"\nProcessing category: {category}")

        # Extract year from category name or histograms
        year = ""
        if "year" in hist_data:
            year = hist_data["year"]
        else:
            year = "13.6 TeV"

        # Process standard (non-neutrino) response variables
        for bin_vars, response_vars, neutrino in zip(
            [BIN_VARIABLES, BIN_VARIABLES_NEUTRINO],
            [RESPONSE_VARIABLES, RESPONSE_VARIABLES_NEUTRINO],
            [False, True],
        ):
            print(
                f"  Processing {'neutrino' if neutrino else 'standard'} response variables..."
            )

            # Build output directory
            output_dir = (
                f"{args.output}/resolution_{category}{'_neutrino' if neutrino else ''}"
            )
            os.makedirs(output_dir, exist_ok=True)

            # Collect histograms for all response variables
            response_hist_dict = {}
            for response_var in response_vars:
                # Map response variable name to histogram name
                hist_key = response_to_hist_name.get(response_var)
                if hist_key and hist_key in hist_data:
                    response_hist_dict[response_var] = hist_data[hist_key]
                else:
                    print(
                        f"  Warning: {response_var} (hist key: {hist_key}) not found in histograms"
                    )

            if not response_hist_dict:
                print(f"  Skipping {category} - no response histograms found")
                continue

            # Get bin variable names from the first response variable histogram
            first_hist = list(response_hist_dict.values())[0]
            bin_var_names = [ax.name for ax in first_hist.axes[:-1]]

            # Plot rho vs pu if available
            rho_hist_key = "h_rho_pu"
            if rho_hist_key in hist_data:
                print(f"  Plotting rho vs pu...")
                plot_rho_vs_pu(hist_data[rho_hist_key], year, category)

            # Compute binned resolution from histograms
            print(f"  Computing binned resolution from histograms...")
            response_tot_dict = compute_binned_resolution_from_hists(
                hist_dict=response_hist_dict,
                bin_var_names=bin_var_names,
                bin_var_configs=bin_vars,
                response_vars=response_vars,
                n_workers=args.workers if args.workers > 1 else 1,
            )

            # Plot resolution vs x variable
            print(f"  Plotting resolution vs x variable...")
            plotters = plot_resolution_vs_x_variable(
                response_types_dict=response_tot_dict,
                bin_var_configs=bin_vars,
                response_vars=response_vars,
                output_dir=output_dir,
                year=year,
            )

            for name, plotter in plotters.items():
                print(f"    Running plotter for bin combination {name}...")
                plotter.run()

    print("\nDone!")


if __name__ == "__main__":
    main()
