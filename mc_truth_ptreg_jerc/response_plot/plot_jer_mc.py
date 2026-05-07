import logging

# Suppress noisy third-party loggers before any imports that trigger them
logging.getLogger("matplotlib").setLevel(logging.WARNING)
logging.getLogger("boost_histogram").setLevel(logging.WARNING)

log = logging.getLogger("plot_jer_mc")

import os
import shutil
import numpy as np
from multiprocessing import Pool
import argparse
import hist
import copy
from scipy.optimize import curve_fit
from scipy.stats import chi2 as chi2_dist
from coffea import util
import json
import yaml


from utils_configs.plot.get_columns_from_files import get_columns_from_files
from utils_configs.plot.HEPPlotter import HEPPlotter

import met_ptreg_performance.helpers as helpers

from mc_truth_ptreg_jerc.response_plot.confidence import Confidence_numpy


def setup_logging(output_dir: str) -> None:
    """Configure root logger with console + rotating file handlers."""
    fmt = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    log.setLevel(logging.DEBUG)

    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    console.setFormatter(fmt)

    log_path = os.path.join(output_dir, "plot_jer_mc.log")
    file_handler = logging.FileHandler(log_path, mode="a")
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(fmt)

    log.addHandler(console)
    log.addHandler(file_handler)
    log.propagate = False
    log.info("Logging to %s", log_path)


def _load_config(config_path):
    """Load and process the YAML configuration file."""
    with open(config_path) as f:
        cfg = yaml.safe_load(f)

    # Convert bin_edges lists to numpy arrays in bin_variable dicts
    for bv_dict_key in ("bin_variables", "bin_variables_neutrino", "bin_variables_mixed"):
        if bv_dict_key not in cfg:
            continue
        for var_cfg in cfg[bv_dict_key].values():
            if "bin_edges" in var_cfg:
                var_cfg["bin_edges"] = np.array(var_cfg["bin_edges"])

    # Apply test mode overrides
    if cfg.get("test", False):
        test_pu = np.array(cfg["test_pu_bins"])
        test_eta = np.array(cfg["test_jet_eta_bins"])
        test_pt = np.array(cfg["test_jet_pt_bins"])
        for bv_dict_key in ("bin_variables", "bin_variables_neutrino", "bin_variables_mixed"):
            if bv_dict_key not in cfg:
                continue
            bv = cfg[bv_dict_key]
            if "Pileup_nPU" in bv:
                bv["Pileup_nPU"]["bin_edges"] = test_pu
            for eta_key in ("MatchedJets_eta", "MatchedJetsNeutrino_eta"):
                if eta_key in bv:
                    bv[eta_key]["bin_edges"] = test_eta
            for pt_key in ("MatchedJets_pt", "MatchedJetsNeutrino_pt"):
                if pt_key in bv:
                    bv[pt_key]["bin_edges"] = test_pt

    # Build bin_variables_mixed as the union of bin_variables and bin_variables_neutrino
    # if it is not explicitly defined in the YAML
    if "bin_variables_mixed" not in cfg:
        cfg["bin_variables_mixed"] = {
            **cfg.get("bin_variables", {}),
            **cfg.get("bin_variables_neutrino", {}),
        }

    # Convert bin_limits lists to tuples in mapping/response variable dicts
    for var_dict_key in (
        "mapping_variables",
        "response_variables",
        "response_variables_neutrino",
        "response_variables_mixed",
    ):
        if var_dict_key not in cfg:
            continue
        for var_cfg in cfg[var_dict_key].values():
            if "bin_limits" in var_cfg:
                var_cfg["bin_limits"] = tuple(var_cfg["bin_limits"])

    # Convert y_lim_resolution list to tuple
    if "y_lim_resolution" in cfg:
        cfg["y_lim_resolution"] = tuple(cfg["y_lim_resolution"])

    return cfg


cfg = None

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
    default=False,
    help="If set, will plot 1d and 2d histograms of the recoil variables",
)
parser.add_argument(
    "-w",
    "--workers",
    type=int,
    default=1,
    help="Number of workers for multiprocessing (default: 1, no multiprocessing)",
)
parser.add_argument(
    "--novars",
    action="store_true",
    help="If true, old save format without saved variations is expected",
    default=False,
)
parser.add_argument(
    "-m",
    "--max-parquet",
    type=int,
    help="Maximum number of parquet files to load",
    default=None,
)
parser.add_argument(
    "-l",
    "--load",
    type=str,
    help="Path to precomputed histograms coffea file to load and plot",
    default=None,
)
parser.add_argument("-o", "--output", type=str, help="Output directory", default="./plot_jer_mc/")
parser.add_argument(
    "-c",
    "--config",
    type=str,
    required=True,
    help="Path to YAML configuration file",
)

args = parser.parse_args()
cfg = _load_config(args.config)


PUPPI_JET_STRING = r"anti-$k_{T}$ R=0.4 (PUPPI)"


def run_plot(plotter):
    """Run a HEPPlotter instance."""
    plotter.run()


def get_color(response_var):
    for key in cfg["plot_settings"]:
        if response_var.endswith(key):
            return cfg["plot_settings"][key]["color"]
    return None


def get_fmt(response_var):
    for key in cfg["plot_settings"]:
        if response_var.endswith(key):
            return cfg["plot_settings"][key]["fmt"]
    return None


def NSC(x, N, S, C, d):
    """
    Relative pT resolution formula.

    Parameters
    ----------
    x : array-like  pT values
    N : float       noise term coefficient
    S : float       stochastic term coefficient
    C : float       constant term
    d : float       additional term coefficient
    """

    return np.sqrt(N * np.abs(N) / (x * x) + S * S * np.power(x, d) + C * C)


# Human-readable model metadata for JSON export and annotations
NSC.formula = "sqrt([0]*abs([0])/(x*x)+[1]*[1]*pow(x,[3])+[2]*[2])"
NSC.param_names = ["[0]", "[1]", "[2]", "[3]"]
NSC.x_name = "x"


def linear_model(x, m, b):
    return m * x + b


linear_model.formula = "[0]*x + [1]"
linear_model.param_names = ["[0]", "[1]"]
linear_model.x_name = "x"


def gaussian_model(x, amplitude, mean, sigma):
    return amplitude * np.exp(-0.5 * ((x - mean) / sigma) ** 2)


gaussian_model.formula = "[0]*exp(-0.5*((x-[1])/[2])^2)"
gaussian_model.param_names = ["[0]", "[1]", "[2]"]
gaussian_model.x_name = "x"


def perform_fit(
    x,
    y,
    y_err,
    fit_function,
    n_params,
    *,
    p0=None,
    bounds=(-np.inf, np.inf),
    absolute_sigma=True,
):
    """
    Perform a (weighted) least-squares fit for generic array data.

    Parameters
    ----------
    x, y : array-like
        Data arrays.
    y_err : array-like or None
        1-sigma uncertainties for y. If None or all non-positive, an unweighted
        fit is performed.
    fit_function : callable
        Function f(x, *params) to fit.
    n_params : int
        Number of fit parameters in fit_function (used for dof calculation).
    p0, bounds, absolute_sigma :
        Passed through to scipy.optimize.curve_fit when applicable.

    Returns
    -------
    dict
        Dictionary containing fit parameters, uncertainties, fit quality metrics,
        and a callable prediction function.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.shape != y.shape:
        raise ValueError("x and y arrays must have the same shape")

    if y_err is None:
        y_err_arr = None
        valid = np.isfinite(x) & np.isfinite(y)
    else:
        y_err_arr = np.asarray(y_err, dtype=float)
        if y_err_arr.shape != y.shape:
            raise ValueError("y_err must have the same shape as y")
        valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(y_err_arr)

    invalid_dict = {
        "fit_function_name": getattr(fit_function, "__name__", str(fit_function)),
        "fit_formula": getattr(fit_function, "formula", None),
        "fit_param_names": getattr(fit_function, "param_names", None),
        "params": np.full(n_params, np.nan),
        "params_err": np.full(n_params, np.nan),
        "cov": None,
        "r_squared": np.nan,
        "chi2": np.nan,
        "dof": 0,
        "p_value": np.nan,
        "n_points": len(x),
        "fit_func": lambda x_val: np.full_like(x_val, np.nan),
        "x_min": np.min(x) if len(x) > 0 else np.nan,
        "x_max": np.max(x) if len(x) > 0 else np.nan,
    }
    if not np.any(valid):
        log.error("No valid points for fitting")
        return invalid_dict

    x = x[valid]
    y = y[valid]
    if y_err_arr is not None:
        y_err_arr = y_err_arr[valid]

    if len(x) < n_params:
        log.error(
            "Not enough valid points for fitting: %d points, %d parameters",
            len(x),
            n_params,
        )
        return invalid_dict

    # If uncertainties are missing or unusable, fall back to unweighted fit.
    use_sigma = (
        y_err_arr is not None and np.any(y_err_arr > 0) and np.all(y_err_arr >= 0)
    )
    sigma = y_err_arr if use_sigma else None

    _curve_fit_kwargs = dict(
        sigma=sigma,
        absolute_sigma=absolute_sigma if use_sigma else False,
        p0=p0,
        bounds=bounds,
        maxfev=100000,
    )
    try:
        popt, pcov = curve_fit(fit_function, x, y, **_curve_fit_kwargs)
    except RuntimeError:
        # Levenberg-Marquardt failed; retry with Trust Region Reflective which
        # handles poorly-scaled problems more robustly.
        try:
            _trf_bounds = bounds if bounds != (-np.inf, np.inf) else (-np.inf, np.inf)
            popt, pcov = curve_fit(
                fit_function,
                x,
                y,
                method="trf",
                **{**_curve_fit_kwargs, "bounds": _trf_bounds},
            )
        except RuntimeError as exc:
            log.error("Fit failed: %s", exc)
            return invalid_dict

    params = np.array(popt, dtype=float)
    params_err = (
        np.sqrt(np.diag(pcov))
        if pcov is not None and pcov.size
        else np.full_like(params, np.nan)
    )

    y_pred = fit_function(x, *params)
    residuals = y - y_pred

    if use_sigma:
        weights = 1.0 / (y_err_arr**2)
        ss_res = np.sum(weights * residuals**2)
        y_mean = np.sum(weights * y) / np.sum(weights)
        ss_tot = np.sum(weights * (y - y_mean) ** 2)
        chi2_stat = ss_res
    else:
        ss_res = np.sum(residuals**2)
        y_mean = np.mean(y)
        ss_tot = np.sum((y - y_mean) ** 2)
        chi2_stat = ss_res

    r_squared = 1.0 - ss_res / ss_tot if ss_tot != 0 else np.nan
    dof = max(len(x) - int(n_params), 1)
    p_value = chi2_dist.sf(chi2_stat, dof)

    log.debug("Fit was successful!")

    return {
        "fit_function_name": getattr(fit_function, "__name__", str(fit_function)),
        "fit_formula": getattr(fit_function, "formula", None),
        "fit_param_names": getattr(fit_function, "param_names", None),
        "params": params,
        "params_err": params_err,
        "cov": pcov,
        "r_squared": r_squared,
        "chi2": chi2_stat,
        "dof": dof,
        "p_value": p_value,
        "n_points": len(x),
        "fit_func": (
            lambda x_val: fit_function(np.asarray(x_val, dtype=float), *params)
        ),
        "x_min": np.min(x) if len(x) > 0 else np.nan,
        "x_max": np.max(x) if len(x) > 0 else np.nan,
    }


def plot_resolution_vs_x_variable(
    response_types_dict,
    response_vars,
    bin_var_configs,
    mapping_dict=None,
    output_dir="",
    year="",
    h_mean_dict=None,
    map_x_variable=False,
    fit_resolution=False,
):
    """
    Plot resolution as a function of the variable marked with resolution_x_variable=True,
    grouped by bins of all other variables (resolution_x_variable=False).
    Overlays multiple response types on the same plot using HEPPlotter.

    Parameters
    ----------
    response_types_dict : dict
        Dictionary mapping response variable names to their resolution_result dicts
        (output from compute_binned_resolution).
        Format: {
            "MatchedJets_ResponseJEC": resolution_result_1,
            "MatchedJets_ResponseRaw": resolution_result_2,
            ...
        }
    response_vars : dict
        Dictionary with configuration for response variables (e.g., RESPONSE_VARIABLES or RESPONSE_VARIABLES_NEUTRINO).
    bin_var_configs : dict
        Configuration dict for binning variables (e.g., BIN_VARIABLES).
    output_dir : str
        Output directory for saving plots (optional).
    year : str
        Year string for annotation.
    h_mean_dict : dict, optional
        Dictionary mapping plot variable names to Mean histograms for mapping x-axis.
    map_x_variable : bool
        If True, use mean values from h_mean_dict for x-axis instead of bin centers.
    mapping_dict : dict, optional
        Dictionary mapping plot variable names to their configurations.
    fit_resolution : bool
        If True, perform a fit to the resolution vs x variable and include fit results in the plot annotation.

    Returns
    -------
    dict
        Dictionary mapping each response variable to its fit results.
        Format: {
            "MatchedJets_ResponseJEC": fit_result_1,
            "MatchedJets_ResponseRaw": fit_result_2,
            ...
        }
    """

    if not response_types_dict:
        raise ValueError("response_types_dict cannot be empty")

    # Pick x_var / y_vars per response, using mapping_dict when mapping the x-axis.
    def _pick_x_var(response_var_name, resolution_result):
        bin_var_names_local = resolution_result["bin_var_names"]

        # Prefer the mapping variable's bin_vars when mapping x.
        if map_x_variable and mapping_dict is not None:
            map_var_name = response_vars[response_var_name].get("map_x_variable")
            if map_var_name and map_var_name in mapping_dict:
                candidate_bin_vars = mapping_dict[map_var_name].get("bin_vars", [])
                x_candidates = [
                    v
                    for v in candidate_bin_vars
                    if bin_var_configs.get(v, {}).get("resolution_x_variable", False)
                ]
                if len(x_candidates) == 1:
                    return x_candidates[0]

        # Fallback: choose from the resolution's own axes
        x_candidates = [
            v
            for v in bin_var_names_local
            if bin_var_configs.get(v, {}).get("resolution_x_variable", False)
        ]
        if len(x_candidates) == 1:
            return x_candidates[0]

        raise ValueError(
            f"Could not uniquely determine x_var for '{response_var_name}'. "
            f"Candidates={x_candidates}, axes={bin_var_names_local}"
        )

    def _name_plot(var_name: str) -> str:
        return bin_var_configs.get(var_name, {}).get("name_plot", var_name)

    # Group response vars by compatible (x_name_plot, y_name_plot tuple) so
    # axes like regular/neutrino that share name_plot get merged onto one plot.
    response_groups = {}
    axis_maps = {}
    for response_var_name, resolution_result in response_types_dict.items():
        bin_var_names_local = resolution_result["bin_var_names"]
        x_var_local = _pick_x_var(response_var_name, resolution_result)

        x_np = _name_plot(x_var_local)
        y_vars_local = [v for v in bin_var_names_local if v != x_var_local]
        y_np = tuple(_name_plot(v) for v in y_vars_local)

        # map merged dims -> concrete axis in this histogram
        by_np = {}
        for v in bin_var_names_local:
            by_np.setdefault(_name_plot(v), []).append(v)

        y_axes_concrete = []
        for dim in y_np:
            candidates = by_np.get(dim, [])
            if len(candidates) != 1:
                raise ValueError(
                    f"Could not uniquely map y-dim '{dim}' for '{response_var_name}'. "
                    f"Candidates={candidates}, axes={bin_var_names_local}"
                )
            y_axes_concrete.append(candidates[0])

        axis_maps[response_var_name] = {
            "x_var": x_var_local,
            "x_np": x_np,
            "y_np": y_np,
            "y_vars": y_axes_concrete,
        }
        response_groups.setdefault((x_np, y_np), []).append(response_var_name)

    all_fit_results = {}

    for (x_np, y_np), group_response_vars in response_groups.items():
        # reference histogram for edges / bin counts in this merged group
        ref_response = group_response_vars[0]
        ref_result = response_types_dict[ref_response]
        ref_axes = axis_maps[ref_response]

        bin_edges_dict = ref_result["bin_edges"]
        bin_var_names_ref = ref_result["bin_var_names"]

        x_var_ref = ref_axes["x_var"]
        y_vars_ref = ref_axes["y_vars"]

        y_vars = [(v, bin_var_names_ref.index(v)) for v in y_vars_ref]

        x_bin_edges = bin_edges_dict[x_var_ref]
        x_bin_centers = (x_bin_edges[:-1] + x_bin_edges[1:]) / 2
        n_x_bins = len(x_bin_centers)

        y_bin_shape = tuple(len(bin_edges_dict[v]) - 1 for v in y_vars_ref)

        plotters = {}

        for y_bin_idx in np.ndindex(y_bin_shape):
            graph_data = {}
            xlabel = None

            for response_var in group_response_vars:
                resolution_result = response_types_dict[response_var]
                resolutions = resolution_result["resolutions"]
                resolutions_uncertainty = resolution_result["resolutions_uncertainty"]
                bin_var_names_local = resolution_result["bin_var_names"]
                axes = axis_maps[response_var]

                x_var_local = axes["x_var"]
                y_vars_local = axes["y_vars"]

                x_var_idx_local = bin_var_names_local.index(x_var_local)
                y_var_idx_local = [bin_var_names_local.index(v) for v in y_vars_local]

                resolutions_for_plot = []
                errors_for_plot = []
                valid_x_indices = []

                for x_idx in range(n_x_bins):
                    full_bin_idx = [0] * len(bin_var_names_local)
                    for i, y_axis_idx in enumerate(y_var_idx_local):
                        full_bin_idx[y_axis_idx] = y_bin_idx[i]
                    full_bin_idx[x_var_idx_local] = x_idx
                    full_bin_idx = tuple(full_bin_idx)

                    if (
                        full_bin_idx in resolutions
                        and resolutions[full_bin_idx] is not None
                    ):
                        resolutions_for_plot.append(resolutions[full_bin_idx])
                        errors_for_plot.append(resolutions_uncertainty[full_bin_idx])
                        valid_x_indices.append(x_idx)
                    else:
                        # put nan for missing points so x-values are still correct if some points are missing
                        resolutions_for_plot.append(np.nan)
                        errors_for_plot.append(np.nan)
                        valid_x_indices.append(x_idx)

                if len(resolutions_for_plot) == 0:
                    continue

                if map_x_variable and h_mean_dict is not None:
                    map_var_name = response_vars[response_var].get("map_x_variable")
                    if not map_var_name or map_var_name not in h_mean_dict:
                        raise ValueError(
                            f"map_x_variable '{map_var_name}' not found in h_mean_dict for response variable '{response_var}'"
                        )
                    h_mean = h_mean_dict[map_var_name]["mean"]
                    # h_mean is expected to share the same axis ordering as the merged dims
                    h_mean_view = h_mean.view()
                    x_values = []
                    for x_idx in valid_x_indices:
                        mean_idx = list(y_bin_idx) + [x_idx]
                        x_values.append(h_mean_view.value[tuple(mean_idx)])
                    x_values = np.asarray(x_values)

                    if xlabel is None:
                        xlabel = (
                            mapping_dict[map_var_name]["label"]
                            if mapping_dict and map_var_name in mapping_dict
                            else map_var_name
                        )
                else:
                    x_values = x_bin_centers[valid_x_indices]
                    if xlabel is None:
                        xlabel = f"{bin_var_configs[x_var_ref]['label']}"

                graph_data[response_var] = {
                    "data": {
                        "x": [x_values, np.zeros_like(x_values)],
                        "y": [
                            np.array(resolutions_for_plot),
                            np.array(errors_for_plot),
                        ],
                    },
                    "style": {
                        "fmt": get_fmt(response_var),
                        "linestyle": "",
                        "color": get_color(response_var),
                        "markersize": 8,
                        "linewidth": 2,
                        "legend_name": response_vars[response_var].get(
                            "legend_name", response_var
                        ),
                        "is_reference": (not map_x_variable)
                        and response_vars[response_var].get("is_reference", False),
                    },
                }

            if not graph_data:
                continue

            if map_x_variable:
                map_var_name = response_vars[response_var].get("map_x_variable")
                name_plot = mapping_dict[map_var_name]["name_plot"]
                filename_parts = ["resolution", f"x_{name_plot}"]
            else:
                filename_parts = ["resolution", f"x_{x_np}"]

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
            output_name = "_".join(filename_parts)

            annotation_text = PUPPI_JET_STRING
            fit_result_string = ""
            for y_var_name, y_idx in y_vars:
                bin_idx = y_bin_idx[y_vars.index((y_var_name, y_idx))]
                low_edge = bin_edges_dict[y_var_name][bin_idx]
                high_edge = bin_edges_dict[y_var_name][bin_idx + 1]
                label = bin_var_configs[y_var_name].get("label", y_var_name)
                annotation_text += f"\n{low_edge} < {label} < {high_edge}"
                fit_result_string += f"_{y_var_name}_{low_edge}to{high_edge}"

            if fit_resolution:
                fit_data = {}
                fit_results = {}
                for response_var, gd in graph_data.items():
                    x_fit = gd["data"]["x"][0]
                    y_fit = gd["data"]["y"][0]
                    y_fit_err = gd["data"]["y"][1]
                    if y_fit_err is not None and np.all(np.asarray(y_fit_err) == 0):
                        y_fit_err = None

                    x_fit = np.asarray(x_fit, dtype=float)
                    pos = np.isfinite(x_fit) & (x_fit > 0)
                    x_fit = x_fit[pos]
                    y_fit = np.asarray(y_fit, dtype=float)[pos]
                    if y_fit_err is not None:
                        y_fit_err = np.asarray(y_fit_err, dtype=float)[pos]

                    if len(x_fit) < len(NSC.param_names):
                        continue

                    fit_res = perform_fit(
                        x_fit,
                        y_fit,
                        y_fit_err,
                        NSC,
                        n_params=len(NSC.param_names),
                        p0=(0.5, 0.5, 0.05, -1.0),
                    )

                    fit_results[f"{response_var}{fit_result_string}"] = fit_res

                    x_fit_fit = np.logspace(
                        np.log10(np.min(x_fit)),
                        np.log10(np.max(x_fit)),
                        100,
                    )
                    y_fit_fit = fit_res["fit_func"](x_fit_fit)
                    fit_data[f"{response_var}_fit"] = {
                        "data": {
                            "x": [x_fit_fit, np.zeros_like(x_fit_fit)],
                            "y": [y_fit_fit, np.zeros_like(y_fit_fit)],
                        },
                        "style": {
                            "linestyle": "--",
                            "fmt": "",
                            "color": get_color(response_var),
                            "linewidth": 2,
                            "appear_in_legend": False,
                        },
                    }

                graph_data.update(fit_data)
                all_fit_results.update(fit_results)

            log.debug("Creating plot for bin combination %s", y_bin_idx)
            plotter = (
                HEPPlotter()
                .set_plot_config(
                    cmstext_loc=2,
                    lumitext=f"{year} (13.6 TeV)",
                    cmstext_font_size=35,
                )
                .set_labels(
                    xlabel=xlabel,
                    ylabel="Jet Energy Resolution",
                    ratio_label="JEC / Resolution",
                )
                .set_data(graph_data, plot_type="graph")
                .set_options(
                    grid=True,
                    legend=True,
                    legend_loc="upper right",
                    x_log=True,
                    split_legend=False,
                    reference_to_den=False,
                    ylim_bottom_value=cfg["y_lim_resolution"][0],
                    ylim_top_value=cfg["y_lim_resolution"][1],
                    set_ylim_ratio=0.5,
                )
                .set_output(f"{output_dir}/{output_name}")
                .add_annotation(
                    0.05,
                    0.7,
                    annotation_text,
                    horizontalalignment="left",
                    verticalalignment="top",
                )
            )

            if fit_resolution:
                for i, (response_var, fit_res) in enumerate(fit_results.items()):
                    chi2_str = f"$\\chi^2$/ndf={fit_res['chi2']:.1f}/{fit_res['dof']}, p={fit_res['p_value']:.2f}"
                    plotter.add_annotation(
                        0.7,
                        0.6 - i * 0.05,
                        chi2_str,
                        color=get_color(response_var.replace(fit_result_string, "")),
                        horizontalalignment="left",
                        verticalalignment="top",
                        fontsize=13,
                    )

            plotters[y_bin_idx] = plotter

        if args.workers == 1:
            for name, plotter in plotters.items():
                log.info(
                    "Plotting resolution vs x variable for bin combination %s...", name
                )
                plotter.run()
        else:
            log.info(
                "Plotting resolution vs x variable for %d bin combinations in parallel with %d workers...",
                len(plotters),
                args.workers,
            )
            with Pool(args.workers) as pool:
                pool.map(run_plot, plotters.values())

    return all_fit_results


def flatten_data(data):
    """
    Flatten all variables in data to a 1D array of length sum(nJets).
    - Jet-level variables (object array of arrays): concatenated directly.
    - Event-level variables (1D flat array): each element repeated nJets[i] times.

    Parameters
    ----------
    data : dict of str -> np.ndarray

    Returns
    -------
    flat_data : dict of str -> np.ndarray
        All variables flattened to the same 1D length.
    """

    def is_ragged(arr):
        return arr.dtype == object and len(arr) > 0 and isinstance(arr[0], np.ndarray)

    # Derive nJets from the first jet-level variable found
    jet_level_key = next(key for key, arr in data.items() if is_ragged(arr))
    nJets = np.array([len(jets) for jets in data[jet_level_key]])

    flat_data = {}
    n_vars = len(data)
    for i, (key, arr) in enumerate(data.items()):
        log.debug("Flattening variable '%s' with shape %s", key, arr.shape)
        if is_ragged(arr):
            flat_data[key] = np.concatenate(arr)
        else:
            flat_data[key] = np.repeat(arr, nJets)
        if (i + 1) % 10 == 0 or (i + 1) == n_vars:
            log.info("Flattening variables: %d/%d done", i + 1, n_vars)

    return flat_data


def create_ND_histo(variables_dict, data, bin_var_configs):
    """
    Create an N-D histogram for each variable in variables_dict, using its bin_vars
    field to determine the binning axes.

    Parameters
    ----------
    variables_dict : dict
        Dictionary mapping variable names to their configurations.
        Each variable must have a "bin_vars" key listing which bin variables to use.
        Optionally, a variable can have a "bin_limits" key: [min, max] for custom axis limits.
    data : dict of str -> np.ndarray
        Dictionary mapping variable names to arrays.
    bin_var_configs : dict
        Configuration dict for binning variables (e.g., BIN_VARIABLES_ALL).

    Returns
    -------
    h_dict : dict
        Maps variable name to hist.Hist object with Count storage.
    h_mean_dict : dict
        Maps variable name to hist.Hist object with Mean storage.
    """
    h_dict = {}
    h_mean_dict = {}
    n_vars = len(variables_dict)

    for i, (var_name, var_cfg) in enumerate(variables_dict.items()):
        log.info("Filling histogram %d/%d: '%s'", i + 1, n_vars, var_name)
        # Get the bin variables for this variable
        bin_var_names = var_cfg.get("bin_vars", [])

        if var_name not in data:
            log.warning("Variable '%s' not found in data, skipping.", var_name)
            continue

        # Check that all required bin variables are available in bin_var_configs
        for bv_name in bin_var_names:
            if bv_name not in bin_var_configs:
                raise ValueError(
                    f"Bin variable '{bv_name}' not found in bin_var_configs"
                )
            if bv_name not in data:
                raise ValueError(f"Bin variable '{bv_name}' not found in data")

        # Build bin axes from the bin_vars
        axes_bin = []
        fill_kwargs_bin = {}
        for bv_name in bin_var_names:
            bv_cfg = bin_var_configs[bv_name]
            axes_bin.append(
                hist.axis.Variable(
                    bv_cfg["bin_edges"], name=bv_name, label=bv_cfg["label"], flow=False
                )
            )
            fill_kwargs_bin[bv_name] = data[bv_name]

        # Build the axis for this variable
        arr = data[var_name]

        # Check if custom bin limits are provided
        if "bin_limits" in var_cfg:
            arr_min, arr_max = var_cfg["bin_limits"]
        else:
            # Auto-compute limits with small padding
            arr_min, arr_max = np.nanmin(arr), np.nanmax(arr)
            padding = (arr_max - arr_min) * 0.001
            arr_min -= padding
            arr_max += padding

        var_axis = hist.axis.Regular(
            var_cfg.get("N_bins", 100),
            arr_min,
            arr_max,
            name=var_name,
            label=var_cfg.get("label", var_name),
            flow=False,
        )

        # Create histograms
        axes_full = axes_bin + [var_axis]
        h = hist.Hist(*axes_full)
        fill_kwargs = copy.deepcopy(fill_kwargs_bin)
        fill_kwargs[var_name] = arr
        h.fill(**fill_kwargs)
        h_dict[var_name] = h

        # Create Mean histogram (only for binning axes)
        if len(axes_bin) > 0:
            h_mean = hist.Hist(*axes_bin, storage=hist.storage.Mean())
            h_mean.fill(**fill_kwargs_bin, sample=arr)
            h_mean_dict[var_name] = h_mean

    return h_dict, h_mean_dict


def compute_means(h_dict, mapping_vars):
    """
    Compute the mean of each mapping variable in bins of its associated bin_vars,
    marginalizing over all other axes.

    Parameters
    ----------
    h_dict : hist.Hist
        The N-D histogram returned by create_ND_histo.
    mapping_vars : dict
        Configuration dict for mapping variables, each with a "bin_vars" key.

    Returns
    -------
    results : dict of str -> dict
        For each plot variable, a dict with:
            "mean"     : hist.Hist with Mean storage, axes = bin_vars of that variable
            "bin_vars" : list of str
    """

    results = {}
    map_vars = list(h_dict.keys())

    for mv in map_vars:
        log.debug("Processing mapping variable: %s", mv)
        map_cfg = mapping_vars[mv]
        active_bin_vars = map_cfg["bin_vars"]

        # Project down to only the active bin axes + this mapping variable axis
        # All other axes are summed over (marginalized)
        h_reduced = h_dict[mv].project(*active_bin_vars, mv)

        # .profile(mv) converts the COUNT histogram into a MEAN histogram
        # by treating the mv axis as the "sample" axis — no numpy needed
        h_mean = h_reduced.profile(mv)

        results[mv] = {
            "mean": h_mean,  # hist.Hist with Mean storage, shape (*active_bin_vars)
            "bin_vars": active_bin_vars,
        }

    return results


def compute_projection(h_dict, mapping_vars, mapping_var_key):
    nd_hist = h_dict[mapping_var_key]
    map_cfg = mapping_vars[mapping_var_key]
    active_bin_vars = map_cfg["bin_vars"] + [mapping_var_key]

    # Project down to only the active bin axes + this mapping variable axis
    # All other axes are summed over (marginalized)
    h_reduced = nd_hist.project(*active_bin_vars)

    return h_reduced


def perform_linear_fit(h_mean_hist):
    """
    Perform a weighted linear fit of a 1D Mean histogram.

    Parameters
    ----------
    h_mean_hist : hist.Hist
        A 1D histogram with Mean storage (e.g. mean of a mapping variable vs its bin variable).

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
    count = np.array(view.count)
    variance = np.array(view.variance)

    if y.shape != x.shape:
        raise ValueError("Histogram x and y arrays must have the same shape")

    y_err = np.where(count > 0, np.sqrt(variance / count), np.inf)
    valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(y_err) & (y_err > 0)
    if not np.any(valid):
        raise ValueError("Histogram contains no valid points for fitting")

    x = x[valid]
    y = y[valid]
    y_err = y_err[valid]

    fit_res = perform_fit(
        x,
        y,
        y_err,
        linear_model,
        n_params=2,
        absolute_sigma=True,
    )

    slope, intercept = fit_res["params"]
    slope_err, intercept_err = fit_res["params_err"]

    return {
        "slope": slope,
        "slope_err": slope_err,
        "intercept": intercept,
        "intercept_err": intercept_err,
        "r_squared": fit_res["r_squared"],
        "chi2": fit_res["chi2"],
        "dof": fit_res["dof"],
        "p_value": fit_res["p_value"],
        "n_points": fit_res["n_points"],
        "fit_func": fit_res["fit_func"],
    }


def _compute_resolution_from_histogram(
    response_var_name, h_response, bin_var_names, response_vars
):
    """
    Compute resolution from an ND histogram for a single response variable.
    Extracts 1D histograms for each bin combination and computes confidence intervals.

    Parameters
    ----------
    response_var_name : str
        Name of the response variable
    h_response : hist.Hist
        ND histogram with axes [bin_var1, bin_var2, ..., response_var]
    bin_var_names : list of str
        Names of the bin variables (in order)
    response_vars : dict
        Configuration dict for response variables

    Returns
    -------
    dict
        Dictionary with keys:
        - "response_var": response variable name
        - "bin_edges": dict mapping bin_var_name to bin edges
        - "bin_var_names": list of bin variable names
        - "resolutions": dict of bin indices -> resolution values
        - "resolution_grid": ndarray with resolution values
    """
    # NOTE: Do not trust the externally provided `bin_var_names` here.
    # Each response variable can have its own per-variable binning axes, and
    # `h_response` encodes the actual axis order/length.
    hist_bin_axes = list(h_response.axes[:-1])
    bin_var_names = [ax.name for ax in hist_bin_axes]

    # Number of bin variable axes
    n_bin_vars = len(hist_bin_axes)
    response_axis_idx = n_bin_vars  # response axis is last

    # Build bin shape
    bin_shape = tuple(len(ax) for ax in hist_bin_axes)

    # Extract bin edges for all bin variables
    bin_edges_dict = {ax.name: np.array(ax.edges) for ax in hist_bin_axes}

    # Get response axis info
    response_axis = h_response.axes[response_axis_idx]
    response_bin_edges = np.array(response_axis.edges)
    response_bin_centers = (response_bin_edges[:-1] + response_bin_edges[1:]) / 2
    response_bin_width = response_bin_edges[1] - response_bin_edges[0]

    resolutions = {}
    resolutions_uncertainty = {}
    resolution_grid = np.full(bin_shape, np.nan)
    gaussian_fits = {}

    log.info("Processing resolution for '%s'...", response_var_name)

    # Get the underlying numpy array from the histogram
    h_view = h_response.view()
    h_variance = h_response.variances()

    # Iterate over all bin combinations of bin variables
    for bin_idx in np.ndindex(bin_shape):
        # Create indexing tuple to extract the 1D histogram along the response axis
        # for this specific bin combination
        slice_tuple = bin_idx + (slice(None),)
        hist_counts = np.array(h_view[slice_tuple])
        hist_variance = np.array(h_variance[slice_tuple])

        failed = False
        # Check if we have enough data (at least 5 events)
        if np.sum(hist_counts) > 5:
            if cfg["gaussian_fit_resolution"]:
                try:
                    x_fit = response_bin_centers.copy()
                    y_fit = hist_counts.astype(float)
                    y_fit_err = np.sqrt(hist_variance)  # Poisson errors

                    if cfg["gaussian_fit_cut_tails"]:
                        # Restrict to the 87% confidence interval around the mean
                        total = y_fit.sum()
                        if total > 0:
                            _, lo_idx, hi_idx = Confidence_numpy(
                                hist=y_fit,
                                bins_mid=x_fit,
                                bin_width=response_bin_width,
                                confLevel=0.87,
                                return_bins=True,
                            )
                            x_fit = x_fit[lo_idx : hi_idx + 1]
                            y_fit = y_fit[lo_idx : hi_idx + 1]
                            y_fit_err = y_fit_err[lo_idx : hi_idx + 1]

                    if len(x_fit) >= 3 and y_fit.sum() > 0:
                        peak_idx = int(np.argmax(y_fit))
                        amplitude_p0 = float(y_fit[peak_idx])
                        mean_p0 = float(x_fit[peak_idx])
                        sigma_p0 = float(response_bin_width * 3)
                        fit_res = perform_fit(
                            x_fit,
                            y_fit,
                            y_err=y_fit_err,
                            fit_function=gaussian_model,
                            n_params=3,
                            p0=(amplitude_p0, mean_p0, sigma_p0),
                            bounds=([0, -np.inf, 0], [np.inf, np.inf, np.inf]),
                        )
                        # Exclude fit_func lambda — lambdas are not picklable and
                        # multiprocessing serializes the return value of this function.
                        gaussian_fits[bin_idx] = {
                            k: v for k, v in fit_res.items() if k != "fit_func"
                        }
                        resolutions[bin_idx] = fit_res["params"][2]
                        resolution_grid[bin_idx] = fit_res["params"][2]
                        resolutions_uncertainty[bin_idx] = fit_res["params_err"][2]
                    else:
                        log.warning(
                            "Not enough valid points for Gaussian fit in bin %s: n_points=%d, total_counts=%.1f",
                            bin_idx,
                            len(x_fit),
                            y_fit.sum(),
                        )
                        failed = True
                except Exception as e:
                    log.error("Failed to fit Gaussian for bin %s: %s", bin_idx, e)
                    failed = True
            else:
                try:
                    resolution = Confidence_numpy(
                        hist=hist_counts,
                        bins_mid=response_bin_centers,
                        bin_width=response_bin_width,
                        confLevel=0.87,
                    )
                    resolutions[bin_idx] = resolution
                    resolution_grid[bin_idx] = resolution
                    resolutions_uncertainty[bin_idx] = (
                        0  # No uncertainty estimate without fit
                    )
                except Exception as e:
                    log.error("Failed to compute resolution for bin %s: %s", bin_idx, e)
                    failed = True
        else:
            log.debug(
                "Skipping bin %s: insufficient data (%d points)",
                bin_idx,
                int(np.sum(hist_counts)),
            )
            failed = True

        if failed:
            log.warning("Setting resolution to NaN for bin %s due to failure", bin_idx)
            resolutions[bin_idx] = np.nan
            resolution_grid[bin_idx] = np.nan
            resolutions_uncertainty[bin_idx] = 0

    return (
        response_var_name,
        {
            "response_var": response_var_name,
            "bin_edges": bin_edges_dict,
            "bin_var_names": bin_var_names,
            "resolutions": resolutions,
            "resolution_grid": resolution_grid,
            "gaussian_fits": gaussian_fits,
            "resolutions_uncertainty": resolutions_uncertainty,
        },
    )


def compute_binned_resolution_from_histograms(h_dict, bin_var_names, response_vars):
    """
    Compute resolution for multiple response variables from pre-computed ND histograms.

    Parameters
    ----------
    h_dict : dict
        Dictionary mapping response_var_name -> hist.Hist (ND histogram)
    bin_var_names : list of str
        Names of the bin variables (in order)
    response_vars : dict
        Configuration dict for response variables, each with a "bin_vars" key

    Returns
    -------
    dict
        Maps response_var_name -> resolution_result dict with keys:
            - "bin_edges": dict mapping bin_var_name to bin edges
            - "bin_var_names": list of bin variable names
            - "resolutions": dict of bin indices -> resolution values
            - "resolution_grid": ndarray with resolution values
    """
    # Prepare arguments for parallel processing
    args_list = [
        (response_var_name, h_response, bin_var_names, response_vars)
        for response_var_name, h_response in h_dict.items()
    ]
    response_tot_dict = {}

    if len(args_list) == 1 or args.workers == 1:
        # Single worker - no need for Pool overhead
        for arg_tuple in args_list:
            log.info(
                "Computing resolution for '%s' without multiprocessing...", arg_tuple[0]
            )
            var_name, result = _compute_resolution_from_histogram(*arg_tuple)
            response_tot_dict[var_name] = result
    else:
        # Multi-worker processing
        log.info(
            "Computing binned resolution for %d response variables in parallel with %d workers...",
            len(args_list),
            args.workers,
        )
        with Pool(args.workers) as pool:
            results = pool.starmap(_compute_resolution_from_histogram, args_list)
            for var_name, result in results:
                response_tot_dict[var_name] = result

    return response_tot_dict


def rebin_histogram(series_dict, target_bins=50, quantile_range=(0.001, 0.999)):
    """
    Rebin all 1D histograms in a dictionary to the same range and fewer bins.
    All histograms will have the same range and binning for proper overlay comparison.

    Parameters
    ----------
    series_dict : dict
        Dictionary mapping variable names to dicts with "data" (hist.Hist) and "style" keys
    target_bins : int
        Approximate number of bins desired after rebinning.
    quantile_range : tuple of float
        Quantile range to trim the axis.

    Returns
    -------
    series_dict : dict
        Dictionary with all histograms rebinned to the same range, preserving style information
    """
    # First, compute combined counts across all histograms to find common range
    all_counts = None
    first_h = None

    for var_name, var_data in series_dict.items():
        h = var_data["data"]
        if h.values().sum() == 0:
            continue

        if first_h is None:
            first_h = h

        if all_counts is None:
            all_counts = np.array(h.values())
        else:
            # Add counts from this histogram
            if len(all_counts) == len(h.values()):
                all_counts = all_counts + np.array(h.values())

    # Check if all histograms are empty
    if all_counts is None or all_counts.sum() == 0:
        log.warning("All histograms are empty, returning original dict")
        return series_dict

    # Get axis from first histogram
    axis = first_h.axes[0]
    centers = axis.centers

    # Compute common range using combined counts
    total = all_counts.sum()
    cdf = np.cumsum(all_counts) / total
    lo_idx = max(0, np.searchsorted(cdf, quantile_range[0]) - 1)
    hi_idx = min(len(centers), np.searchsorted(cdf, quantile_range[1]) + 1)
    n_bins_in_range = hi_idx - lo_idx

    # Compute rebin factor for the common range.
    # Use floor division so the output has *at least* target_bins bins.
    # round() could overshoot (e.g. round(75/50)=2 → 38 bins instead of 50).
    rebin_factor = max(1, n_bins_in_range // target_bins)
    remainder = n_bins_in_range % rebin_factor
    if remainder != 0:
        hi_idx = min(len(centers), hi_idx + (rebin_factor - remainder))

    lo_edge = axis.edges[lo_idx]
    hi_edge = axis.edges[hi_idx]

    # Apply the same rebinning to all histograms
    rebinned_dict = {}
    for var_name, var_data in series_dict.items():
        h = var_data["data"]
        h_rebinned = h[lo_edge * 1j : hi_edge * 1j : rebin_factor * 1j]

        rebinned_dict[var_name] = {
            "data": h_rebinned,
            "style": var_data["style"],
        }

    return rebinned_dict


def plot_variable_slices(
    h_dict,
    variables_dict,
    bin_var_configs,
    output_dir="",
    category="",
    year="",
    var_type="mapping",
    gaussian_fit_results=None,
):
    """
    Plot all 1D histogram slices from ND histograms.
    Overlays all variables on the same plot. Works for both plot and response variables.
    Variables are grouped by their bin structure (different bin_vars are plotted separately).

    Parameters
    ----------
    h_dict : dict
        Dictionary mapping variable names to ND hist.Hist objects
    variables_dict : dict
        Configuration dict for variables (plot or response)
    bin_var_configs : dict
        Configuration dict for bin variables
    output_dir : str
        Output directory for plots
    category : str
        Category name for output filenames
    year : str
        Year string for annotation
    var_type : str
        Type of variables: "mapping" (default) or "response"
    gaussian_fit_results : dict, optional
        Maps variable name -> dict of bin_idx -> fit_result (from _compute_resolution_from_histogram).
        When provided and GAUSSIAN_FIT_RESOLUTION is True, the Gaussian fit curve is overlaid on
        each histogram slice.

    Returns
    -------
    None (saves plots to disk)
    """
    type_label = (
        f"{var_type} variables"
        if var_type in ("mapping", "plot")
        else "response variables"
    )
    log.info("Plotting %s slices for %s...", type_label, category)

    # Group variables by their bin structure (they may have different bin_vars)
    variables_by_bin_axes = {}
    for var_name, h_var in h_dict.items():
        n_bin_axes = h_var.ndim - 1
        bin_var_names_for_var = tuple([h_var.axes[i].name for i in range(n_bin_axes)])

        if bin_var_names_for_var not in variables_by_bin_axes:
            variables_by_bin_axes[bin_var_names_for_var] = {}
        variables_by_bin_axes[bin_var_names_for_var][var_name] = h_var

    plotters = []

    # Process each bin structure group separately
    for bin_var_names_tuple, group_dict in variables_by_bin_axes.items():
        bin_var_names = list(bin_var_names_tuple)
        n_bin_axes = len(bin_var_names)

        if n_bin_axes == 0:
            # No binning variables, just plot all variables in this group together
            series_dict = {}
            for var_name, h_var in group_dict.items():
                var_cfg = variables_dict[var_name]
                series_dict[var_name] = {
                    "data": h_var,
                    "style": {
                        "color": get_color(var_name),
                        "legend_name": var_cfg.get("legend_name", var_name),
                    },
                }

            # Create output filename
            output_name = f"{output_dir}/histo_{variables_dict[var_name]['name_plot']}_slice_{category}"

            log.debug(
                "Creating plot for unbinned variables: %s...", list(group_dict.keys())
            )
            # Plot using HEPPlotter
            plotter = (
                HEPPlotter()
                .set_plot_config(lumitext=f"{year} (13.6 TeV)" if year else "")
                .set_data(series_dict, plot_type="1d")
                .set_labels(
                    xlabel=f"{list(group_dict.values())[0].axes[-1].label}",
                    ylabel="Events",
                )
                .set_output(output_name)
                .set_options(grid=True, legend=True)
            )
            plotters.append(plotter)
            continue

        # Get the first histogram in this group to determine bin shape
        first_h = list(group_dict.values())[0]
        bin_shape = tuple(len(first_h.axes[i]) for i in range(n_bin_axes))

        # Iterate over all bin combinations
        for bin_idx in np.ndindex(bin_shape):
            # Build annotation text from bin ranges
            annotation_text = f"{PUPPI_JET_STRING}\n"
            for i, bin_i in enumerate(bin_idx):
                bin_var_name = bin_var_names[i]
                bin_cfg = bin_var_configs[bin_var_name]
                bin_edges = bin_cfg["bin_edges"]
                low_edge = bin_edges[bin_i]
                high_edge = bin_edges[bin_i + 1]
                label = bin_cfg.get("label", bin_var_name)
                if isinstance(low_edge, np.int64) and isinstance(high_edge, np.int64):
                    annotation_text += f"{low_edge} <= {label} < {high_edge}\n"
                else:
                    annotation_text += f"{low_edge:.2f} < {label} < {high_edge:.2f}\n"

            # Collect histograms for this bin combination
            series_dict = {}
            for var_name, h_var in group_dict.items():
                var_cfg = variables_dict[var_name]

                # Extract 1D histogram for this bin combination
                slice_tuple = bin_idx + (slice(None),)
                h_1d = h_var[slice_tuple]

                series_dict[var_name] = {
                    "data": h_1d,
                    "style": {
                        "color": get_color(var_name),
                        "legend_name": var_cfg.get("legend_name", var_name),
                        "linewidth": 2,
                    },
                }

            # check if all histograms in this group have empty data for this bin combination
            if all(
                np.sum(h_var["data"].values()) == 0 for h_var in series_dict.values()
            ):
                log.warning(
                    "All histograms in group %s are empty for bin combination %s, skipping plot.",
                    bin_var_names,
                    bin_idx,
                )
                continue

            if variables_dict[var_name].get("rebin_for_plotting", False):
                series_dict = rebin_histogram(series_dict, 50, (0.01, 0.99))

            # Create output filename with bin ranges
            filename_parts = [
                f"histo_{variables_dict[var_name]['name_plot']}_slice",
                category,
            ]
            for i, bin_i in enumerate(bin_idx):
                bin_var_name = bin_var_names[i]
                bin_cfg = bin_var_configs[bin_var_name]
                bin_edges = bin_cfg["bin_edges"]
                low_edge_str = f"{bin_edges[bin_i]}".replace(".", "p").replace("-", "m")
                high_edge_str = f"{bin_edges[bin_i + 1]}".replace(".", "p").replace(
                    "-", "m"
                )
                filename_parts.append(
                    f"{bin_cfg['name_plot']}_{low_edge_str}to{high_edge_str}"
                )

            output_name = f"{output_dir}/{'_'.join(filename_parts)}"

            # Get the axis label from the first variable in this group
            axis_label = list(group_dict.values())[0].axes[-1].label

            log.debug(
                "Creating plot for bin combination %s with variables: %s...",
                bin_idx,
                list(group_dict.keys()),
            )
            # Plot using HEPPlotter
            plotter = (
                HEPPlotter()
                .set_plot_config(
                    lumitext=f"{year} (13.6 TeV)" if year else "", figsize=(13, 13)
                )
                .set_data(series_dict, plot_type="1d")
                .set_labels(xlabel=axis_label, ylabel="Events")
                .set_output(output_name)
                .set_options(
                    grid=True,
                    legend=True,
                    y_log=True,
                    split_legend=False,
                    ylim_bottom_value=1,
                )
                .add_annotation(
                    0.05,
                    0.95,
                    annotation_text,
                    coord_type="axes",
                    verticalalignment="top",
                    fontsize=20,
                )
            )

            # Overlay Gaussian fit curves when fit results are available
            if gaussian_fit_results is not None and cfg["gaussian_fit_resolution"]:
                for i, var_name in enumerate(group_dict):

                    var_fits = gaussian_fit_results.get(var_name, {})
                    fit_res = var_fits.get(bin_idx)
                    if fit_res is None or np.any(np.isnan(fit_res["params"])):
                        continue
                    h_1d = series_dict[var_name]["data"]
                    axis = h_1d.axes[0]
                    x_fine = np.linspace(axis.edges[0], axis.edges[-1], 300)
                    amplitude = fit_res["params"][0]
                    mu = fit_res["params"][1]
                    sigma = fit_res["params"][2]

                    y_fine = gaussian_model(x_fine, amplitude, mu, sigma)
                    fit_label = (
                        f"$\\sigma={abs(sigma):.3f} \\pm {fit_res['params_err'][2]:.3f}$\n"
                        f"$\\chi^2/\\mathrm{{ndf}}={fit_res['chi2']:.0f}/{fit_res['dof']}$, "
                        f"$p={fit_res['p_value']:.2f}$"
                    )
                    plotter.add_curve(
                        x_fine,
                        y_fine,
                        color=get_color(var_name),
                        linestyle="--",
                        linewidth=2,
                    )
                    plotter.add_annotation(
                        0.7,
                        0.75 - 0.1 * i,
                        fit_label,
                        coord_type="axes",
                        verticalalignment="top",
                        fontsize=20,
                        color=get_color(var_name),
                    )

                    if cfg["gaussian_fit_cut_tails"]:
                        # plot the lines corresponding to the fit range
                        plotter.add_line(
                            "v", x=fit_res["x_min"], color=get_color(var_name), linestyle="dotted"
                        )
                        plotter.add_line(
                            "v", x=fit_res["x_max"], color=get_color(var_name), linestyle="dotted"
                        )

            plotters.append(plotter)

    # Run all plotters
    if args.workers == 1:
        for plotter in plotters:
            log.debug("Plotting %s slice...", type_label)
            plotter.run()
    else:
        log.info(
            "Plotting %s slices in parallel with %d workers...",
            type_label,
            args.workers,
        )
        with Pool(args.workers) as pool:
            pool.map(run_plot, plotters)


def profile_means(h_mean_dict, mapping_vars):
    """
    Compute the mean of each mapping variable in bins of its associated bin_vars,
    marginalizing over all other axes, using hist's built-in Mean storage.

    Parameters
    ----------
    h_mean_dict : dict of str -> hist.Hist
        Histograms with Mean storage for each mapping variable, returned by create_ND_histo.
    mapping_vars : dict
        Configuration dict for mapping variables, each with a "bin_vars" key.

    Returns
    -------
    results : dict of str -> dict
        For each mapping variable, a dict with:
            "mean"     : hist.Hist with Mean storage, axes = bin_vars of that variable
            "bin_vars" : list of str
    """

    results = {}

    for mv, h_mean in h_mean_dict.items():
        log.info("Processing mapping variable: %s", mv)
        map_cfg = mapping_vars[mv]
        active_bin_vars = map_cfg["bin_vars"]

        # Project down to only the active bin axes
        # The Mean storage already contains the mean values and uncertainties
        h_reduced = h_mean.project(*active_bin_vars)

        # set zeros to nan
        h_reduced.values()[h_reduced.values() == 0] = np.nan

        results[mv] = {
            "mean": h_reduced,  # hist.Hist with Mean storage, shape (*active_bin_vars)
            "bin_vars": active_bin_vars,
        }

    return results


def save_plotted_data(output_path, data_dict):
    """
    Save plotted data (histograms, results, etc.) to a coffea file.

    Parameters
    ----------
    output_path : str
        Path to save the coffea file
    data_dict : dict
        Dictionary containing all plotted data
    """

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    util.save(data_dict, output_path)
    log.info("Saved plotted data to %s", output_path)


def load_plotted_data(input_path):
    """
    Load plotted data (histograms, results, etc.) from a coffea file.

    Parameters
    ----------
    input_path : str
        Path to the coffea file

    Returns
    -------
    dict
        Dictionary containing all plotted data
    """

    data_dict = util.load(input_path)
    log.info("Loaded plotted data from %s", input_path)
    return data_dict


def plot_mapping_variable_linear_fit(
    h_dict,
    results,
    mapping_vars,
    year,
    category,
    output_dir,
    *,
    var_name,
    bin_var_configs,
):
    var_cfg = mapping_vars[var_name]
    x_var_name = var_cfg["bin_vars"][0]
    x_cfg = bin_var_configs[x_var_name]

    x_label = x_cfg["label"]
    x_name = x_cfg["name_plot"]
    y_label = var_cfg["label"]
    y_name = var_cfg["name_plot"]
    y_mean_label = f"$< {y_label.strip('$')} >$"

    # 2D histogram of mapping variable vs its bin variable
    var_2d = compute_projection(h_dict, mapping_vars, var_name)
    (
        HEPPlotter()
        .set_plot_config(lumitext=f"{year} (13.6 TeV)")
        .set_options(legend=False, cbar_log=False)
        .set_output(f"{output_dir}/{y_name}_vs_{x_name}_2d_{category}")
        .set_data({"data": {"data": var_2d, "style": {}}}, plot_type="2d")
        .set_labels(x_label, y_label)
    ).run()

    fit_results = perform_linear_fit(results[var_name]["mean"])
    log.debug("Linear fit results: %s", fit_results)

    x_axis = results[var_name]["mean"].axes[0]
    x_lin = np.linspace(x_axis.edges[0], x_axis.edges[-1], 100)
    mean_dict = {
        "data": {
            "data": {
                "x": [(x_axis.edges[:-1] + x_axis.edges[1:]) / 2, x_axis.widths / 2],
                "y": [
                    results[var_name]["mean"].view().value,
                    np.sqrt(
                        results[var_name]["mean"].view().variance
                        / np.where(
                            results[var_name]["mean"].view().count > 0,
                            results[var_name]["mean"].view().count,
                            np.inf,
                        )
                    ),
                ],
            },
            "style": {"fmt": "o"},
        },
        "linear fit": {
            "data": {
                "x": [x_lin, np.zeros(100)],
                "y": [fit_results["fit_func"](x_lin), np.zeros(100)],
            },
            "style": {"linestyle": "-", "fmt": ""},
        },
    }
    (
        HEPPlotter()
        .set_plot_config(lumitext=f"{year} (13.6 TeV)")
        .set_options(set_ylim=False)
        .set_output(f"{output_dir}/{y_name}_vs_{x_name}_fit_{category}")
        .set_data(mean_dict, plot_type="graph")
        .set_labels(x_label, y_mean_label)
        .add_annotation(
            x=0.05,
            y=0.7,
            s=f"$\\chi^2/ndf=${fit_results['chi2']:.0f}/{fit_results['dof']}\n"
            + f"$p$-value={fit_results['p_value']:.3f}\n",
        )
    ).run()

    return fit_results


def save_fit_results(fit_results, output_path):
    """
    Save fit results to a json file.
    """
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)

    def _to_builtin(obj):
        if obj is None:
            return None
        if isinstance(obj, (str, int, float, bool)):
            return obj
        if isinstance(obj, (list, tuple)):
            return [_to_builtin(x) for x in obj]
        if isinstance(obj, dict):
            return {str(k): _to_builtin(v) for k, v in obj.items()}
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, (np.floating, np.integer)):
            return obj.item()
        return str(obj)

    serializable = {}
    for response_var, fit_result in (fit_results or {}).items():
        if fit_result is None:
            continue

        # Drop the callable from JSON; keep metadata and numerical results.
        fr = dict(fit_result)
        fr.pop("fit_func", None)

        serializable[response_var] = _to_builtin(fr)

    with open(output_path, "w") as f:
        json.dump(serializable, f, indent=4)
    log.info("Saved fit results to %s", output_path)


def save_txt_resolution(
    fit_results,
    response_var_name,
    response_var_cfg,
    bin_var_configs,
    linear_fit_maps,
    mapping_vars,
    output_dir,
    year,
):
    """
    Save NSC resolution fit results in JERC text format.

    Header format:  { N_bin_vars  BinVar1  BinVar2  1  XVar  formula  Resolution}
    Row format:     bin1_lo  bin1_hi  ...  (2+n_params)  x_min  x_max  p0  p1  ...

    Bin variables with "txt_map_to" are mapped through the corresponding linear fit
    (e.g. Pileup_nPU mu edges → Rho edges). Non-mapped variables appear first.
    """

    def get_txt_name(var_name):
        bin_var_cfg = bin_var_configs[var_name]
        if "txt_name" in bin_var_cfg:
            return bin_var_cfg["txt_name"]
        mapped_to = bin_var_cfg.get("txt_map_to")
        if mapped_to and mapped_to in mapping_vars:
            return mapping_vars[mapped_to].get("txt_name", mapped_to)
        return var_name

    x_var = next(
        (k for k, v in bin_var_configs.items() if v.get("resolution_x_variable")),
        None,
    )

    if x_var is None:
        log.warning(
            "No resolution_x_variable found for %s, skipping txt.", response_var_name
        )
        return

    prefix = f"{response_var_name}_"
    matching = {k: v for k, v in fit_results.items() if k.startswith(prefix)}
    if not matching:
        log.warning("No fit results for %s, skipping txt output.", response_var_name)
        return

    # Detect which bin variables are actually in the keys. In "mixed" mode,
    # bin_var_configs contains both regular and neutrino variables, but a given
    # response var only uses a subset of them.
    first_key_suffix = next(iter(matching.keys()))[len(response_var_name) :]
    vars_in_key = [
        k for k in bin_var_configs if k != x_var and f"_{k}_" in first_key_suffix
    ]
    non_mapped = [k for k in vars_in_key if "txt_name" in bin_var_configs[k]]
    mapped = [k for k in vars_in_key if "txt_map_to" in bin_var_configs[k]]
    y_vars_txt = non_mapped + mapped
    y_vars_hist = vars_in_key

    first_res = next(iter(matching.values()))
    formula = first_res.get("fit_formula", "unknown")
    n_fit_params = len(first_res["params"])
    n_vals = 2 + n_fit_params

    n_bin_vars = len(y_vars_txt)
    header = (
        "{ "
        + f"{n_bin_vars} "
        + " ".join(get_txt_name(v) for v in y_vars_txt)
        + f" 1 {get_txt_name(x_var)} "
        + formula
        + " Resolution}"
    )

    rows = []
    for key, fit_res in matching.items():
        suffix = key[len(response_var_name) :]
        bin_edges_raw = {}
        parse_ok = True
        for y_var in y_vars_hist:
            marker = f"_{y_var}_"
            idx = suffix.find(marker)
            if idx == -1:
                parse_ok = False
                break
            val_start = idx + len(marker)
            val_end = len(suffix)
            for other_var in y_vars_hist:
                other_idx = suffix.find(f"_{other_var}_", val_start)
                if other_idx != -1 and other_idx < val_end:
                    val_end = other_idx
            try:
                lo_str, hi_str = suffix[val_start:val_end].split("to", 1)
                bin_edges_raw[y_var] = (float(lo_str), float(hi_str))
            except ValueError:
                parse_ok = False
                break

        if not parse_ok:
            log.warning("Could not parse bin edges from key '%s', skipping row.", key)
            continue

        row_cols = []
        build_ok = True
        for y_var in y_vars_txt:
            if y_var not in bin_edges_raw:
                build_ok = False
                break
            lo, hi = bin_edges_raw[y_var]
            bin_var_cfg = bin_var_configs[y_var]
            if "txt_map_to" in bin_var_cfg:
                fit_map = (linear_fit_maps or {}).get(bin_var_cfg["txt_map_to"])
                if fit_map is not None:
                    lo, hi = fit_map["fit_func"](np.array([lo, hi]))
            row_cols.extend([lo, hi])

        if not build_ok:
            log.warning("Missing bin variable in row for key '%s', skipping.", key)
            continue

        rows.append(
            row_cols
            + [n_vals, fit_res["x_min"], fit_res["x_max"]]
            + list(fit_res["params"])
        )

    if not rows:
        log.warning("No valid rows for %s, skipping txt output.", response_var_name)
        return

    rows.sort(key=lambda r: r[: 2 * n_bin_vars])

    def _fmt(val):
        if isinstance(val, (int, np.integer)):
            return str(int(val))
        return f"{float(val):g}"

    year_tag = cfg["year_map"].get(year, year)
    jet_type = response_var_cfg.get("txt_jet_name", "AK4PFPuppi")
    filename = f"Run3{year_tag}_V1_NSC_MC_PtResolution_{jet_type}.txt"
    output_path = os.path.join(output_dir, filename)
    os.makedirs(output_dir, exist_ok=True)

    with open(output_path, "w") as f:
        f.write(header + "\n")
        for row in rows:
            f.write("  ".join(_fmt(v) for v in row) + "\n")
    log.info("Saved resolution txt to %s", output_path)


def plot_mapping_variable_histograms(*, category, year, h_dict):
    if args.histo or cfg["histograms_map"]:
        plot_var_output_dir = f"{args.output}/histograms_mapping_variables_{category}"
        os.makedirs(plot_var_output_dir, exist_ok=True)
        plot_variable_slices(
            h_dict=h_dict,
            variables_dict=cfg["mapping_variables"],
            bin_var_configs=cfg["bin_variables_mixed"],
            output_dir=plot_var_output_dir,
            category=category,
            year=year,
            var_type="mapping",
        )


def get_bin_vars_for_mode(mode):
    if mode == "regular":
        return cfg["bin_variables"]
    elif mode == "neutrino":
        return cfg["bin_variables_neutrino"]
    elif mode == "mixed":
        return cfg["bin_variables_mixed"]
    else:
        raise ValueError(f"Unknown mode: {mode}")


def process_response_type(
    *,
    category,
    year,
    mode,
    response_h_dict,
    response_tot_dict,
    available_response_vars,
    bin_vars,
    results_for_mapping,
    linear_fit_maps=None,
):
    suffix = f"_{mode}"

    if args.histo or cfg["histograms_response"]:
        response_histogram_dir = (
            f"{args.output}/histograms_resolution_{category}{suffix}"
        )
        os.makedirs(response_histogram_dir, exist_ok=True)
        # Extract Gaussian fit results from the precomputed resolution dicts
        if cfg["gaussian_fit_resolution"] and response_tot_dict:
            gaussian_fit_results = {
                var_name: result.get("gaussian_fits", {})
                for var_name, result in response_tot_dict.items()
            }
            save_fit_results(
                gaussian_fit_results,
                f"{args.output}/gaussian_fit_results_{category}{suffix}.json",
            )

        else:
            gaussian_fit_results = None

        plot_variable_slices(
            h_dict=response_h_dict,
            variables_dict=available_response_vars,
            bin_var_configs=bin_vars,
            output_dir=response_histogram_dir,
            category=category,
            year=year,
            var_type="response",
            gaussian_fit_results=gaussian_fit_results,
        )

    if cfg["resolution_vs_pt_gen"]:
        # Plot resolution vs x variable for response variables (bin centers)
        output_dir = f"{args.output}/resolution_ptgen_{category}{suffix}"
        os.makedirs(output_dir, exist_ok=True)
        plot_resolution_vs_x_variable(
            response_types_dict=response_tot_dict,
            bin_var_configs=bin_vars,
            response_vars=available_response_vars,
            mapping_dict=cfg["mapping_variables"],
            output_dir=output_dir,
            year=year,
        )

    if cfg["resolution_vs_pt_reco"]:
        # Plot resolution vs mapped x variable (means) + fit + save fit results
        output_dir = f"{args.output}/resolution_ptreco_{category}{suffix}"
        os.makedirs(output_dir, exist_ok=True)
        fit_results = plot_resolution_vs_x_variable(
            response_types_dict=response_tot_dict,
            bin_var_configs=bin_vars,
            response_vars=available_response_vars,
            mapping_dict=cfg["mapping_variables"],
            output_dir=output_dir,
            year=year,
            h_mean_dict=results_for_mapping,
            map_x_variable=True,
            fit_resolution=True,
        )
        save_fit_results(
            fit_results,
            f"{args.output}/resolution_fit_results_{category}{suffix}.json",
        )

        for resp_var_name, resp_var_cfg in available_response_vars.items():
            save_txt_resolution(
                fit_results=fit_results,
                response_var_name=resp_var_name,
                response_var_cfg=resp_var_cfg,
                bin_var_configs=bin_vars,
                linear_fit_maps=linear_fit_maps,
                mapping_vars=cfg["mapping_variables"],
                output_dir=args.output,
                year=year,
            )


def main():
    # make output dir if it doesn't exist
    os.makedirs(args.output, exist_ok=True)
    shutil.copy2(args.config, os.path.join(args.output, os.path.basename(args.config)))
    setup_logging(args.output)

    # If loading from precomputed file, skip data loading and computation
    if args.load:
        log.info("Loading precomputed data from %s...", args.load)
        loaded_data = load_plotted_data(args.load)
        year = loaded_data["year"]
        category = loaded_data["category"]
        cat_data = loaded_data["histogram_data"]
        h_dict = cat_data["h_dict"]
        results = cat_data["results"]

        # Load linear fit results from file if available, otherwise recompute
        if "linear_fit_maps" in cat_data and False:
            linear_fit_maps = cat_data["linear_fit_maps"]
            mapped_bin_edges = cat_data.get("mapped_bin_edges", {})
        else:
            linear_fit_maps = {}
            mapped_bin_edges = {}
            for lf_var_name, lf_var_cfg in cfg["mapping_variables"].items():
                if not lf_var_cfg.get("linear_fit", False):
                    continue
                linear_fit_maps[lf_var_name] = plot_mapping_variable_linear_fit(
                    h_dict=h_dict,
                    results=results,
                    mapping_vars=cfg["mapping_variables"],
                    year=year,
                    category=category,
                    output_dir=args.output,
                    var_name=lf_var_name,
                    bin_var_configs=cfg["bin_variables_mixed"],
                )
                mapped_bin_edges[lf_var_name] = linear_fit_maps[lf_var_name][
                    "fit_func"
                ](cfg["bin_variables"][lf_var_cfg["bin_vars"][0]]["bin_edges"])

        plot_mapping_variable_histograms(category=category, year=year, h_dict=h_dict)

        # Process response variables from loaded data
        for mode in ["regular", "neutrino", "mixed"]:
            if mode not in cat_data["response_data"]:
                continue

            resp_data = cat_data["response_data"][mode]
            response_h_dict = resp_data["response_h_dict"]
            response_tot_dict = resp_data["response_tot_dict"]
            available_response_vars = resp_data["available_response_vars"]

            bin_vars = get_bin_vars_for_mode(mode)

            process_response_type(
                category=category,
                year=year,
                mode=mode,
                response_h_dict=response_h_dict,
                response_tot_dict=response_tot_dict,
                available_response_vars=available_response_vars,
                bin_vars=bin_vars,
                results_for_mapping=results,
                linear_fit_maps=linear_fit_maps,
            )
        return

    inputfiles_data = [
        os.path.join(args.input_dir, file)
        for file in os.listdir(args.input_dir)
        if file.endswith(".coffea")
    ]

    log.info("Loading the columns...")
    cat_col, total_datasets_list = get_columns_from_files(
        inputfiles_data,
        "nominal",
        None,
        True,
        args.novars,
        max_num_parquet_files=args.max_parquet,
    )

    # get year from dataset
    year = ""
    for dataset in total_datasets_list:
        # dataset are of the shape name_year_yearsuffix or name_year
        y = helpers.extract_year_tag(dataset)

        year = ", ".join([year, y]) if year else y

    log.info("Total datasets found: %s", total_datasets_list)
    log.info("Year: %s", year)

    for category, col_var in cat_col.items():
        log.info("Processing category '%s': flattening data...", category)
        col_var_flatten = flatten_data(col_var)

        # Process plot variables (contains both regular and neutrino versions)
        log.info("Building mapping variable histograms for category '%s'...", category)
        h_dict, h_mean_dict = create_ND_histo(
            variables_dict=cfg["mapping_variables"],
            data=col_var_flatten,
            bin_var_configs=cfg["bin_variables_mixed"],
        )

        results = profile_means(h_mean_dict, cfg["mapping_variables"])

        # Fit mapping variables that have a linear_fit defined
        linear_fit_maps = {}
        mapped_bin_edges = {}
        for lf_var_name, lf_var_cfg in cfg["mapping_variables"].items():
            if not lf_var_cfg.get("linear_fit", False):
                continue
            linear_fit_maps[lf_var_name] = plot_mapping_variable_linear_fit(
                h_dict=h_dict,
                results=results,
                mapping_vars=cfg["mapping_variables"],
                year=year,
                category=category,
                output_dir=args.output,
                var_name=lf_var_name,
                bin_var_configs=cfg["bin_variables_mixed"],
            )
            mapped_bin_edges[lf_var_name] = linear_fit_maps[lf_var_name]["fit_func"](
                cfg["bin_variables"][lf_var_cfg["bin_vars"][0]]["bin_edges"]
            )

        # Store category data for saving
        category_data = {
            "h_dict": h_dict,
            "h_mean_dict": h_mean_dict,
            "results": results,
            "linear_fit_maps": linear_fit_maps,
            "mapped_bin_edges": mapped_bin_edges,
            "response_data": {},
        }

        # Process response variables separately for neutrino and non-neutrino
        for response_vars, mode in (
            zip(
                [cfg["response_variables"], cfg["response_variables_neutrino"]],
                ["regular", "neutrino"],
            )
            if not cfg["mixed_mode"]
            else zip(
                [cfg["response_variables_mixed"]],
                ["mixed"],
            )
        ):
            # Determine which bin_vars to use based on neutrino flag
            bin_vars = get_bin_vars_for_mode(mode)

            bin_var_names = list(bin_vars.keys())

            # Filter response variables to only those available in data
            available_response_vars = {
                response_var: response_vars[response_var]
                for response_var in response_vars
                if response_var in col_var_flatten
            }

            if not available_response_vars:
                log.warning("No response variables found for mode=%s, skipping.", mode)
                continue

            # Create ND histograms for response variables
            log.info(
                "Building response histograms for mode='%s', category='%s'...",
                mode,
                category,
            )
            response_h_dict, response_h_mean_dict = create_ND_histo(
                variables_dict=available_response_vars,
                data=col_var_flatten,
                bin_var_configs=bin_vars,
            )

            # Compute resolutions from the ND histograms
            response_tot_dict = compute_binned_resolution_from_histograms(
                h_dict=response_h_dict,
                bin_var_names=bin_var_names,
                response_vars=available_response_vars,
            )

            # Store response data for this category
            response_type_key = mode
            category_data["response_data"][response_type_key] = {
                "response_h_dict": response_h_dict,
                "response_tot_dict": response_tot_dict,
                "available_response_vars": available_response_vars,
                "bin_vars": bin_vars,
            }

        # Save all plotted data to coffea file
        output_coffea = os.path.join(args.output, f"plotted_data_{category}.coffea")
        save_plotted_data(
            output_coffea,
            {
                "year": year,
                "category": category,
                "histogram_data": category_data,
            },
        )

        # Plot
        plot_mapping_variable_histograms(category=category, year=year, h_dict=h_dict)
        for response_type, resp_data in category_data["response_data"].items():

            process_response_type(
                category=category,
                year=year,
                mode=mode,
                response_h_dict=resp_data["response_h_dict"],
                response_tot_dict=resp_data["response_tot_dict"],
                available_response_vars=resp_data["available_response_vars"],
                bin_vars=resp_data["bin_vars"],
                results_for_mapping=results,
                linear_fit_maps=linear_fit_maps,
            )

    log.info("All done!")


if __name__ == "__main__":

    main()
