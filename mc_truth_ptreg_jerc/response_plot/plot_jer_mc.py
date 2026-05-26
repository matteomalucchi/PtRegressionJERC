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


def _load_config(config_path, test_override=False):
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
    if test_override:
        test_pu = np.array(cfg["test_pu_bins"]) if "test_pu_bins" in cfg else None
        test_eta = np.array(cfg["test_jet_eta_bins"]) if "test_jet_eta_bins" in cfg else None
        test_pt = np.array(cfg["test_jet_pt_bins"]) if "test_jet_pt_bins" in cfg else None
        for bv_dict_key in ("bin_variables", "bin_variables_neutrino", "bin_variables_mixed"):
            if bv_dict_key not in cfg:
                continue
            bv = cfg[bv_dict_key]
            for key in bv:
                key_lower = key.lower()
                if "npu" in key_lower and test_pu is not None:
                    bv[key]["bin_edges"] = test_pu
                elif "eta" in key_lower and test_eta is not None:
                    bv[key]["bin_edges"] = test_eta
                elif "pt" in key_lower and test_pt is not None:
                    bv[key]["bin_edges"] = test_pt

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
parser.add_argument(
    "--refit",
    action="store_true",
    default=False,
    help=(
        "When used with --load, merge histogram bins using the bin edges from the current YAML "
        "config (if coarser than what is stored) and recompute the Gaussian fit before plotting. "
        "Bins in the coffea file that are a strict sub-group of the config edges are merged; "
        "if the config requests finer binning than available a warning is emitted and the axis "
        "is left unchanged."
    ),
)
parser.add_argument(
    "--test",
    action="store_true",
    default=False,
    help="Run in test mode with reduced bin arrays (overrides the 'test' flag in the YAML config)",
)

args = parser.parse_args()
cfg = _load_config(args.config, test_override=args.test)


PUPPI_JET_STRING = r"anti-$k_{T}$ R=0.4 (PUPPI)"

# When True: removes grid, puts CMS text inside plots (cmstext_loc=2),
# removes log scale from histograms, and removes resolution annotations from histograms.
DP_NOTE_PLOTS = True


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

        y_bin_label_dict = build_bin_label_dict(y_bin_shape, y_vars_ref, bin_edges_dict)

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

                eta_max = response_vars[response_var].get("eta_max")
                if eta_max is not None:
                    bin_edges_local = resolution_result["bin_edges"]
                    if any(
                        "eta" in vn.lower()
                        and min(
                            abs(bin_edges_local[vn][y_bin_idx[i]]),
                            abs(bin_edges_local[vn][y_bin_idx[i] + 1]),
                        ) >= eta_max
                        for i, vn in enumerate(y_vars_local)
                    ):
                        continue

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

                if all(np.isnan(v) for v in resolutions_for_plot):
                    log.warning(
                        "All resolutions are NaN for %s in bin %s, skipping plot",
                        response_var,
                        y_bin_label_dict[y_bin_idx]["label"],
                    )
                    continue

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

            annotation_text = PUPPI_JET_STRING
            fit_result_string = ""
            for y_var_name, _ in y_vars:
                low_edge, high_edge = y_bin_label_dict[y_bin_idx]["edges"][y_var_name]
                low_edge_str = (f"{low_edge}" if isinstance(low_edge, (int, np.integer)) else f"{low_edge}".replace(".", "p")).replace("-", "m")
                high_edge_str = (f"{high_edge}" if isinstance(high_edge, (int, np.integer)) else f"{high_edge}".replace(".", "p")).replace("-", "m")
                filename_parts.append(
                    f"{bin_var_configs[y_var_name]['name_plot']}_{low_edge_str}to{high_edge_str}"
                )
                label = bin_var_configs[y_var_name].get("label", y_var_name)
                annotation_text += "\n" + _format_bin_annotation_line(low_edge, high_edge, label)
                fit_result_string += f"_{y_var_name}_{low_edge}to{high_edge}"
                
            output_name = "_".join(filename_parts)

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

            log.debug("Creating plot for bin combination %s", y_bin_label_dict[y_bin_idx]["label"])

            ylim_bottom = cfg["y_lim_resolution"][0]
            ylim_top = cfg["y_lim_resolution"][1]
            has_ratio = any(
                gd.get("style", {}).get("is_reference", False)
                for gd in graph_data.values()
            )
            data_y_max = ylim_bottom
            for gd in graph_data.values():
                y_vals = np.asarray(gd["data"]["y"][0], dtype=float)
                y_err = gd["data"]["y"][1]
                if y_err is not None:
                    y_err = np.asarray(y_err, dtype=float)
                    y_vals = y_vals + np.abs(y_err)
                finite = y_vals[np.isfinite(y_vals)]
                if len(finite):
                    data_y_max = max(data_y_max, float(finite.max()))
            ylim_top_base = cfg["y_lim_resolution"][1]
            if has_ratio:
                # Extra headroom for the ratio panel; annotation moves to 0.65
                ylim_top = max(ylim_top_base + 0.2, data_y_max + 0.10)
                annotation_y = 0.67
            else:
                if data_y_max > ylim_top * 0.6:
                    ylim_top = ylim_top_base + 0.2  
                annotation_y = 0.75

            plotter = (
                HEPPlotter()
                .set_plot_config(
                    cmstext="Simulation\nPreliminary",
                    cmstext_loc=2,
                    lumitext=f"{year} (13.6 TeV)",
                    cmstext_font_size=30,
                )
                .set_labels(
                    xlabel=xlabel,
                    ylabel="Jet Energy Resolution",
                    ratio_label="JEC / Resolution",
                )
                .set_data(graph_data, plot_type="graph")
                .set_options(
                    grid=not DP_NOTE_PLOTS,
                    legend=True,
                    legend_loc="upper right",
                    x_log=True,
                    split_legend=False,
                    reference_to_den=False,
                    ylim_bottom_value=ylim_bottom,
                    ylim_top_value=ylim_top,
                    set_ylim_ratio=0.5,
                )
                .set_output(f"{output_dir}/{output_name}")
                .add_annotation(
                    0.05,
                    annotation_y,
                    annotation_text,
                    horizontalalignment="left",
                    verticalalignment="top",
                )
            )

            if fit_resolution:
                for i, (response_var, fit_res) in enumerate(fit_results.items()):
                    chi2_str = f"$\\chi^2$/ndf={fit_res['chi2']:.0f}/{fit_res['dof']}, p={fit_res['p_value']:.2f}"
                    plotter.add_annotation(
                        0.69,
                        0.65 - i * 0.07,
                        chi2_str,
                        color=get_color(response_var.replace(fit_result_string, "")),
                        horizontalalignment="left",
                        verticalalignment="top",
                        fontsize=15,
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
    h_mean_bin_var_dict : dict
        Maps bin variable name to a 1D hist.Hist with Mean storage of that
        variable in its own bins (used to compute mean bin-var value per bin).
    """
    h_dict = {}
    h_mean_dict = {}
    h_mean_bin_var_dict = {}
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

            # For each bin variable, build a 1D Mean histogram of that variable
            # in its own bins so callers can use mean bin values instead of centers.
            for bv_name, bv_axis in zip(bin_var_names, axes_bin):
                if bv_name not in h_mean_bin_var_dict:
                    h_mean_bv = hist.Hist(bv_axis, storage=hist.storage.Mean())
                    h_mean_bv.fill(
                        **{bv_name: fill_kwargs_bin[bv_name]},
                        sample=fill_kwargs_bin[bv_name],
                    )
                    h_mean_bin_var_dict[bv_name] = h_mean_bv

    return h_dict, h_mean_dict, h_mean_bin_var_dict


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
        if mv not in mapping_vars:
            log.warning(
                "Mapping variable '%s' found in histogram dict but not in mapping_vars config; skipping.",
                mv,
            )
            continue
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
    existing_axes = {ax.name for ax in nd_hist.axes}
    active_bin_vars = [v for v in map_cfg["bin_vars"] + [mapping_var_key] if v in existing_axes]

    # Project down to only the active bin axes + this mapping variable axis
    # All other axes are summed over (marginalized)
    h_reduced = nd_hist.project(*active_bin_vars)

    return h_reduced


def perform_linear_fit(h_mean_hist, x_values=None):
    """
    Perform a weighted linear fit of a 1D Mean histogram.

    Parameters
    ----------
    h_mean_hist : hist.Hist
        A 1D histogram with Mean storage (e.g. mean of a mapping variable vs its bin variable).
    x_values : array-like, optional
        If provided, use these values as the x coordinates instead of the bin centers.
        Must have the same length as the number of bins. NaN entries are excluded from the fit.

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

    x = np.array(x_values) if x_values is not None else np.array(x_axis.centers)
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


def merge_nd_histogram_bins(h, bin_var_configs_new):
    """
    Merge histogram bins along any axis whose name appears in *bin_var_configs_new* and
    whose new bin edges are a coarser (proper sub-set) version of the axis's current edges.

    Works for both ``hist.storage.Count`` and ``hist.storage.Mean`` storages.
    Axes whose names are absent from *bin_var_configs_new* (e.g. the response-variable
    axis) are left unchanged.

    Parameters
    ----------
    h : hist.Hist
        Histogram to rebin.  May have any number of axes.
    bin_var_configs_new : dict
        Mapping from variable name to config dict that must contain a ``bin_edges`` key
        with the desired (possibly coarser) bin edges.

    Returns
    -------
    h_new : hist.Hist
        New histogram with merged bins.
    any_rebinned : bool
        ``True`` if at least one axis was actually merged.
    """
    # Check via the view dtype: Mean storage returns a structured array with named fields
    # ('value', 'count'); Count storage returns a plain float array.
    is_mean = h.view().dtype.names is not None

    arr = h.view().copy()
    new_axes = []
    any_rebinned = False

    for axis_idx, ax in enumerate(h.axes):
        var_name = ax.name

        if var_name not in bin_var_configs_new:
            new_axes.append(ax)
            continue

        orig_edges = np.array(ax.edges)
        new_edges = np.array(bin_var_configs_new[var_name]["bin_edges"])

        if len(new_edges) >= len(orig_edges):
            log.warning(
                "Variable '%s': new config has %d edges but original has %d — "
                "cannot increase binning, axis left unchanged.",
                var_name,
                len(new_edges),
                len(orig_edges),
            )
            new_axes.append(ax)
            continue

        # Verify every new edge is present in the original edges (within tolerance).
        all_match = all(
            np.any(np.isclose(ne, orig_edges, rtol=1e-6, atol=1e-10)) for ne in new_edges
        )
        if not all_match:
            log.warning(
                "New bin edges for '%s' are not a subset of the original edges — "
                "axis left unchanged.\n  original: %s\n  requested: %s",
                var_name,
                orig_edges,
                new_edges,
            )
            new_axes.append(ax)
            continue

        # Map each new lower edge to its index in orig_edges.
        group_starts = [
            int(np.argmin(np.abs(orig_edges - ne))) for ne in new_edges[:-1]
        ]

        if is_mean:
            # Mean storage view: structured array whose dtype must be preserved
            # exactly (boost_histogram requires all fields, e.g.
            # ('count', 'value', '_sum_of_deltas_squared')).
            # _sum_of_deltas_squared is summed additively — this underestimates the
            # true within-group variance by the between-bin term, but keeps it
            # non-zero so that downstream fits (which require y_err > 0) still work.
            log.warning(
                "Merging bins for variable '%s' with Mean storage: this underestimates the true variance within merged bins, but keeps it non-zero for downstream fits.",
                var_name,
            )
            orig_dtype = arr.dtype
            counts = arr["count"].copy()
            sum_vals = counts * arr["value"]

            merged_counts = np.add.reduceat(counts, group_starts, axis=axis_idx)
            merged_sum_vals = np.add.reduceat(sum_vals, group_starts, axis=axis_idx)

            new_arr = np.zeros(merged_counts.shape, dtype=orig_dtype)
            new_arr["count"] = merged_counts
            new_arr["value"] = np.where(
                merged_counts > 0, merged_sum_vals / merged_counts, 0.0
            )
            if "_sum_of_deltas_squared" in orig_dtype.names:
                new_arr["_sum_of_deltas_squared"] = np.add.reduceat(
                    arr["_sum_of_deltas_squared"].copy(), group_starts, axis=axis_idx
                )
            arr = new_arr
        else:
            arr = np.add.reduceat(arr, group_starts, axis=axis_idx)

        new_axes.append(
            hist.axis.Variable(new_edges, name=ax.name, label=ax.label, flow=False)
        )
        any_rebinned = True
        log.info(
            "Merged bins for '%s': %d → %d bins",
            var_name,
            len(orig_edges) - 1,
            len(new_edges) - 1,
        )

    if not any_rebinned:
        return h, False

    storage = hist.storage.Mean() if is_mean else hist.storage.Double()
    h_new = hist.Hist(*new_axes, storage=storage)
    if is_mean:
        # Assign field-by-field to handle dtype mismatches between coffea files
        # saved with different boost_histogram versions (e.g. 2-field vs 3-field
        # Mean storage). Fields absent in arr are left at zero.
        view = h_new.view()
        for field in arr.dtype.names:
            if field in view.dtype.names:
                view[field] = arr[field]
    else:
        h_new.view()[:] = arr
    return h_new, True


def squeeze_single_bin_axes(h, bin_var_configs):
    """Remove histogram axes that are bin variables with exactly one bin.

    Only axes whose names appear in *bin_var_configs* are candidates so that
    the response/mapping variable axis is never touched.

    Returns the squeezed histogram and a list of removed axis names.
    """
    slice_list = []
    removed = []
    for ax in h.axes:
        if ax.name in bin_var_configs and len(ax) == 1:
            slice_list.append(0)
            removed.append(ax.name)
        else:
            slice_list.append(slice(None))
    if not removed:
        return h, []
    if len(removed) == h.ndim:
        log.warning(
            "Skipping squeeze: all axes are single-bin (%s); would produce a scalar, not a histogram.",
            removed,
        )
        return h, []
    log.info("Squeezing single-bin axes %s from histogram", removed)
    return h[tuple(slice_list)], removed


def _bin_excluded_by_eta_max(bin_idx, bin_var_names, bin_edges_dict, eta_max):
    """Return True if this bin lies entirely outside |eta| < eta_max for any eta axis."""
    if eta_max is None:
        return False
    for i, vn in enumerate(bin_var_names):
        if "eta" in vn.lower():
            edges = bin_edges_dict[vn]
            low_abs = abs(edges[bin_idx[i]])
            high_abs = abs(edges[bin_idx[i] + 1])
            if min(low_abs, high_abs) >= eta_max:
                return True
    return False


def build_bin_label_dict(bin_shape, bin_var_names, bin_edges_dict):
    """Build a dict mapping bin_idx tuples to a dict with a human-readable label
    string and per-variable (low, high) edge pairs.

    Each entry: {"label": str, "edges": {var_name: (low, high)}}
    """
    result = {}
    for bin_idx in np.ndindex(bin_shape):
        parts = []
        edges_per_var = {}
        for i, var_name in enumerate(bin_var_names):
            edges = bin_edges_dict[var_name]
            low, high = edges[bin_idx[i]], edges[bin_idx[i] + 1]
            edges_per_var[var_name] = (low, high)
            def _is_int_like(v):
                return isinstance(v, (int, np.integer)) or (isinstance(v, float) and v.is_integer())
            if _is_int_like(low) and _is_int_like(high):
                parts.append(f"{var_name}: [{int(low)}, {int(high)})")
            else:
                parts.append(f"{var_name}: [{low:.2f}, {high:.2f})")
        result[bin_idx] = {"label": ", ".join(parts), "edges": edges_per_var}
    return result


def _format_bin_annotation_line(low_edge, high_edge, label):
    """Format one bin variable as a single annotation line (no newline)."""
    if low_edge.is_integer() and high_edge.is_integer():
        return f"{int(low_edge)} <= {label} < {int(high_edge)}"
    return f"{low_edge:.2f} < {label} < {high_edge:.2f}"


def _parse_additional_uncertainty(value, response_var_name=""):
    """Parse the per-response-variable ``additional_uncertainty`` config option.

    The option may be provided as a number or as a string (per request). Empty,
    None, or unparseable values default to 0.
    """
    if value is None or value == "":
        return 0.0
    try:
        return float(value)
    except (TypeError, ValueError):
        log.warning(
            "Could not parse additional_uncertainty=%r for '%s' as float; using 0",
            value,
            response_var_name,
        )
        return 0.0


def _apply_response_normalization_and_extra_uncertainty(
    response_var_name,
    resolutions,
    resolutions_uncertainty,
    resolution_grid,
    response_means,
    response_means_uncertainty,
    response_vars,
    gaussian_fits=None,
):
    """Apply optional per-bin transformations to the computed resolution values.

    All transformations are configured per response variable in the YAML
    config. When applied, the order matches the conventional one used by
    upstream JER tooling:

    1. ``normalize_by_response_mean`` (bool, default False): divide each
       resolution by the mean of the response (sigma / <R>), with the
       propagated uncertainty
           eerel = erel * sqrt((sigma_err/sigma)^2 + (mean_err/mean)^2)
    2. ``add_min_err`` (float or string, default 0.0): add a minimal
       uncertainty floor in quadrature
           new_err = sqrt(err^2 + add_min_err^2)
       (e.g. 0.0005 for a 0.05% floor).
    3. ``add_chi2_ndof`` (bool, default False): multiply the uncertainty by
       ``sqrt(chi2 / ndf)`` whenever the per-bin Gaussian fit's reduced
       chi-square exceeds 1, leaving it unchanged otherwise. Requires the
       Gaussian fit to be enabled and successful for that bin.
    4. ``additional_uncertainty`` (float or string, default 0.0): a relative
       uncertainty added in quadrature
           new_err = sqrt(err^2 + (additional_uncertainty * resolution)^2)
    """
    var_cfg = response_vars.get(response_var_name, {})
    normalize = bool(var_cfg.get("normalize_by_response_mean", False))
    add_unc = _parse_additional_uncertainty(
        var_cfg.get("additional_uncertainty", 0.0), response_var_name
    )
    add_min_err = _parse_additional_uncertainty(
        var_cfg.get("add_min_err", 0.0), response_var_name
    )
    add_chi2_ndof = bool(var_cfg.get("add_chi2_ndof", False))
    gaussian_fits = gaussian_fits or {}

    if not normalize and add_unc <= 0 and add_min_err <= 0 and not add_chi2_ndof:
        return

    if normalize:
        log.info(
            "Normalizing resolution by mean of response for '%s'",
            response_var_name,
        )
    if add_min_err > 0:
        log.info(
            "Adding minimal uncertainty floor %.4f in quadrature for '%s'",
            add_min_err,
            response_var_name,
        )
    if add_chi2_ndof:
        log.info(
            "Scaling uncertainty by sqrt(chi2/ndf) where chi2/ndf > 1 for '%s'",
            response_var_name,
        )
    if add_unc > 0:
        log.info(
            "Adding extra relative uncertainty %.4f in quadrature for '%s'",
            add_unc,
            response_var_name,
        )

    for bin_idx in list(resolutions.keys()):
        sigma = resolutions[bin_idx]
        sigma_err = resolutions_uncertainty.get(bin_idx, 0.0)
        if sigma is None or np.isnan(sigma):
            continue

        if normalize:
            mean = response_means.get(bin_idx, np.nan)
            mean_err = response_means_uncertainty.get(bin_idx, 0.0)
            if (
                mean is None
                or np.isnan(mean)
                or mean == 0
                or sigma / mean < 0
            ):
                resolutions[bin_idx] = np.nan
                resolution_grid[bin_idx] = np.nan
                resolutions_uncertainty[bin_idx] = 0.0
                continue
            sigma_rel = sigma / mean
            if sigma > 0:
                sigma_rel_err = abs(sigma_rel) * np.sqrt(
                    (sigma_err / sigma) ** 2 + (mean_err / mean) ** 2
                )
            else:
                sigma_rel_err = abs(sigma_rel * mean_err / mean) if mean != 0 else 0.0
            sigma = sigma_rel
            sigma_err = sigma_rel_err

        if add_min_err > 0:
            sigma_err = float(np.sqrt(sigma_err**2 + add_min_err**2))

        if add_chi2_ndof:
            fit_res = gaussian_fits.get(bin_idx)
            if fit_res is not None:
                chi2 = fit_res.get("chi2", np.nan)
                dof = fit_res.get("dof", 0) or 0
                if dof > 0 and np.isfinite(chi2):
                    chi2_per_ndof = float(chi2) / float(dof)
                    if chi2_per_ndof > 1:
                        sigma_err = float(sigma_err * np.sqrt(chi2_per_ndof))

        if add_unc > 0:
            sigma_err = float(np.sqrt(sigma_err**2 + (add_unc * sigma) ** 2))

        resolutions[bin_idx] = sigma
        resolution_grid[bin_idx] = sigma
        resolutions_uncertainty[bin_idx] = sigma_err


def _apply_response_mean_method(
    *,
    response_var_name,
    h_response,
    h_mean,
    response_means,
    response_means_uncertainty,
    gaussian_fits,
    bin_shape,
    response_axis,
    response_vars,
):
    """Optionally override the per-bin response mean using a configurable method.

    Controlled by the per-response-variable config key ``response_mean_method``,
    which may be one of:

    - ``"gaussian_fit"``: use the mean from the per-bin Gaussian fit
      (``params[1]`` / ``params_err[1]``). Bins where the Gaussian fit is
      unavailable or invalid get NaN.
    - ``"histogram"``: compute the first moment from the binned 1D slice of
      the ND response histogram (``sum(y*x) / sum(y)``), with the standard
      error of the mean as the uncertainty.
    - ``"mean_storage"``: read the mean and its uncertainty directly from a
      ``hist.storage.Mean()`` histogram filled with the un-binned response
      values (as produced by ``create_ND_histo``). The uncertainty is the
      standard error of the mean ``sqrt(variance / count)``.

    When the option is unset (or set to ``"auto"`` / None) the response means
    that were stored inline during ``_compute_resolution_from_histogram`` are
    left untouched (Gaussian fit mean when the resolution came from a
    successful Gaussian fit, binned moments otherwise).
    """
    var_cfg = response_vars.get(response_var_name, {})
    method = var_cfg.get("response_mean_method")
    if method in (None, "", "auto"):
        return

    method_str = str(method).strip().lower()
    valid_methods = ("gaussian_fit", "histogram", "mean_storage")
    if method_str not in valid_methods:
        log.warning(
            "Unknown response_mean_method=%r for '%s'; valid: %s. "
            "Keeping default response means.",
            method,
            response_var_name,
            valid_methods,
        )
        return

    log.info(
        "Computing response mean for '%s' using method '%s'",
        response_var_name,
        method_str,
    )

    if method_str == "gaussian_fit":
        for bin_idx in np.ndindex(bin_shape):
            fit_res = gaussian_fits.get(bin_idx)
            if fit_res is None or "params" not in fit_res:
                response_means[bin_idx] = np.nan
                response_means_uncertainty[bin_idx] = 0.0
                continue
            params = np.asarray(fit_res["params"], dtype=float)
            params_err = np.asarray(fit_res.get("params_err", []), dtype=float)
            if params.size < 2 or np.any(np.isnan(params[:2])):
                response_means[bin_idx] = np.nan
                response_means_uncertainty[bin_idx] = 0.0
                continue
            response_means[bin_idx] = float(params[1])
            response_means_uncertainty[bin_idx] = (
                float(params_err[1]) if params_err.size > 1 and np.isfinite(params_err[1]) else 0.0
            )
        return

    if method_str == "mean_storage":
        if h_mean is None:
            log.warning(
                "response_mean_method='mean_storage' requested for '%s' but no "
                "Mean-storage histogram was provided; keeping default response means.",
                response_var_name,
            )
            return
        view = h_mean.view()
        try:
            value_arr = np.asarray(view.value)
            count_arr = np.asarray(view.count)
            variance_arr = np.asarray(view.variance)
        except AttributeError:
            log.warning(
                "Provided Mean-storage histogram for '%s' does not expose "
                "value/count/variance; keeping default response means.",
                response_var_name,
            )
            return
        if value_arr.shape != bin_shape:
            log.warning(
                "Mean-storage histogram shape %s does not match response bin "
                "shape %s for '%s'; keeping default response means.",
                value_arr.shape,
                bin_shape,
                response_var_name,
            )
            return
        with np.errstate(invalid="ignore", divide="ignore"):
            mean_err_arr = np.where(
                (count_arr > 0) & (variance_arr > 0),
                np.sqrt(variance_arr / np.where(count_arr > 0, count_arr, 1)),
                0.0,
            )
        for bin_idx in np.ndindex(bin_shape):
            if count_arr[bin_idx] > 0:
                response_means[bin_idx] = float(value_arr[bin_idx])
                response_means_uncertainty[bin_idx] = float(mean_err_arr[bin_idx])
            else:
                response_means[bin_idx] = np.nan
                response_means_uncertainty[bin_idx] = 0.0
        return

    # method_str == "histogram": vectorised binned moments over the response axis
    h_view = np.asarray(h_response.view())
    edges = np.asarray(response_axis.edges)
    centers = (edges[:-1] + edges[1:]) / 2
    bc_shape = (1,) * (h_view.ndim - 1) + (-1,)
    bc = centers.reshape(bc_shape)

    with np.errstate(invalid="ignore", divide="ignore"):
        ntot = np.sum(h_view, axis=-1)
        safe_ntot = np.where(ntot > 0, ntot, 1)
        means = np.where(ntot > 0, np.sum(h_view * bc, axis=-1) / safe_ntot, np.nan)
        diff_sq = (bc - means[..., None]) ** 2
        mu2 = np.where(
            ntot > 0,
            np.sum(h_view * diff_sq, axis=-1) / safe_ntot,
            np.nan,
        )
        mean_err = np.where(
            (ntot > 0) & (mu2 > 0),
            np.sqrt(np.where(mu2 > 0, mu2, 0) / safe_ntot),
            0.0,
        )

    for bin_idx in np.ndindex(bin_shape):
        m = means[bin_idx]
        if np.isnan(m):
            response_means[bin_idx] = np.nan
            response_means_uncertainty[bin_idx] = 0.0
        else:
            response_means[bin_idx] = float(m)
            response_means_uncertainty[bin_idx] = float(mean_err[bin_idx])


def _compute_resolution_from_histogram(
    response_var_name, h_response, bin_var_names, response_vars, h_mean=None
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
        Configuration dict for response variables. Each entry may include:
        - "normalize_by_response_mean" (bool, default False): when True, divide
          the per-bin resolution by the mean of the response (sigma / <R>) and
          propagate the uncertainty accordingly.
        - "additional_uncertainty" (float or string, default 0.0): a relative
          uncertainty added in quadrature to the resolution uncertainty as
          new_err = sqrt(err^2 + (additional_uncertainty * resolution)^2).
        - "add_min_err" (float or string, default 0.0): a minimal uncertainty
          floor added in quadrature as new_err = sqrt(err^2 + add_min_err^2).
        - "add_chi2_ndof" (bool, default False): when True, multiply the
          uncertainty by sqrt(chi2/ndf) wherever the per-bin Gaussian fit's
          reduced chi-square is greater than 1.
        - "response_mean_method" (str, default None / "auto"): how to compute
          the per-bin mean of the response. One of:
            * "gaussian_fit"  - mean from the per-bin Gaussian fit
            * "histogram"     - first moment of the binned 1D response slice
            * "mean_storage"  - mean from a hist.storage.Mean() histogram
              filled with the un-binned response values (provided via
              ``h_mean``).
          When unset, the default behaviour stores the Gaussian fit mean when
          the resolution itself comes from a successful Gaussian fit and the
          binned histogram mean otherwise.
    h_mean : hist.Hist or None, optional
        Optional Mean-storage histogram (axes = bin variables only) holding
        the mean of the response per bin. Used by
        ``response_mean_method="mean_storage"`` and ignored otherwise.

    Returns
    -------
    dict
        Dictionary with keys:
        - "response_var": response variable name
        - "bin_edges": dict mapping bin_var_name to bin edges
        - "bin_var_names": list of bin variable names
        - "resolutions": dict of bin indices -> resolution values (possibly
          normalized by the response mean)
        - "resolution_grid": ndarray with resolution values
        - "resolutions_uncertainty": dict of bin indices -> uncertainty on
          resolution (possibly normalized and with extra uncertainty added)
        - "response_means": dict of bin indices -> mean of response per bin
        - "response_means_uncertainty": dict of bin indices -> uncertainty on
          the mean of response per bin
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

    bin_label_dict = build_bin_label_dict(bin_shape, bin_var_names, bin_edges_dict)

    # Get response axis info
    response_axis = h_response.axes[response_axis_idx]
    response_bin_edges = np.array(response_axis.edges)
    response_bin_centers = (response_bin_edges[:-1] + response_bin_edges[1:]) / 2
    response_bin_width = response_bin_edges[1] - response_bin_edges[0]

    resolutions = {}
    resolutions_uncertainty = {}
    resolution_grid = np.full(bin_shape, np.nan)
    gaussian_fits = {}
    rebin_infos = {}
    ci_intervals = {}
    response_means = {}
    response_means_uncertainty = {}

    log.info("Processing resolution for '%s'...", response_var_name)

    eta_max = response_vars[response_var_name].get("eta_max")

    # Get the underlying numpy array from the histogram
    h_view = h_response.view()
    h_variance = h_response.variances()

    # Iterate over all bin combinations of bin variables
    for bin_idx in np.ndindex(bin_shape):
        if _bin_excluded_by_eta_max(bin_idx, bin_var_names, bin_edges_dict, eta_max):
            log.debug(
                "Skipping bin %s for '%s': outside eta_max=%.2f",
                bin_label_dict[bin_idx]["label"], response_var_name, eta_max,
            )
            resolutions[bin_idx] = np.nan
            resolution_grid[bin_idx] = np.nan
            resolutions_uncertainty[bin_idx] = 0
            response_means[bin_idx] = np.nan
            response_means_uncertainty[bin_idx] = 0
            continue

        # Create indexing tuple to extract the 1D histogram along the response axis
        # for this specific bin combination
        slice_tuple = bin_idx + (slice(None),)
        hist_counts = np.array(h_view[slice_tuple])
        hist_variance = np.array(h_variance[slice_tuple])

        failed = False
        rebin_info_for_idx = None
        # Check if we have enough data
        if np.sum(hist_counts) > cfg.get("min_events_for_fit", 50):
            x_fit = response_bin_centers.copy()
            y_fit = hist_counts.astype(float)
            y_fit_err = np.sqrt(hist_variance)  # Poisson errors
            bin_width_for_fit = response_bin_width

            if response_vars[response_var_name].get("rebin_for_plotting", False):
                # Build a temporary 1D histogram so rebin_histogram can operate on it
                h_1d_slice = hist.Hist(
                    hist.axis.Variable(
                        response_bin_edges,
                        name=response_axis.name,
                        label=response_axis.label,
                        flow=False,
                    )
                )
                h_1d_slice.view()[:] = hist_counts
                h_1d_rebinned, rebin_info = rebin_histogram(h_1d_slice)
                if rebin_info is not None:
                    x_fit = h_1d_rebinned.axes[0].centers
                    y_fit = h_1d_rebinned.values().astype(float)
                    y_fit_err = np.sqrt(h_1d_rebinned.variances())
                    rebin_info_for_idx = rebin_info
                    bin_width_for_fit = response_bin_width * rebin_info[2]

            if rebin_info_for_idx is not None:
                rebin_infos[bin_idx] = rebin_info_for_idx

            if cfg["gaussian_fit_resolution"]:
                try:
                    _cut_cfg = cfg.get("gaussian_fit_cut_tails", False)
                    if isinstance(_cut_cfg, list):
                        ci_values_to_try = [float(v) if v not in (None, False) else None for v in _cut_cfg]
                    elif _cut_cfg is False or _cut_cfg is None:
                        ci_values_to_try = [None]
                    else:
                        ci_values_to_try = [float(_cut_cfg)]

                    max_rel_err = cfg.get("gaussian_fit_max_sigma_rel_err", 1.0)
                    x_fit_base = x_fit.copy()
                    y_fit_base = y_fit.copy()
                    y_fit_err_base = y_fit_err.copy()
                    fit_res = None
                    fit_accepted = False

                    # Collect fit results for all CI values, then select the best by p-value.
                    fit_candidates = []  # (ci, fit_res, x_fit_try, y_fit_try, y_fit_err_try, is_valid)
                    for ci in ci_values_to_try:
                        x_fit_try = x_fit_base.copy()
                        y_fit_try = y_fit_base.copy()
                        y_fit_err_try = y_fit_err_base.copy()

                        if ci is not None:
                            total = y_fit_try.sum()
                            if total > 0:
                                _, lo_idx, hi_idx = Confidence_numpy(
                                    hist=y_fit_try,
                                    bins_mid=x_fit_try,
                                    bin_width=bin_width_for_fit,
                                    confLevel=ci,
                                    return_bins=True,
                                )
                                x_fit_try = x_fit_try[lo_idx : hi_idx + 1]
                                y_fit_try = y_fit_try[lo_idx : hi_idx + 1]
                                y_fit_err_try = y_fit_err_try[lo_idx : hi_idx + 1]

                        if len(x_fit_try) < 3 or y_fit_try.sum() <= 0:
                            log.debug(
                                "CI=%s: not enough points in bin %s, skipping",
                                ci, bin_label_dict[bin_idx]["label"],
                            )
                            continue

                        peak_idx = int(np.argmax(y_fit_try))
                        fit_res_try = perform_fit(
                            x_fit_try,
                            y_fit_try,
                            y_err=y_fit_err_try,
                            fit_function=gaussian_model,
                            n_params=3,
                            p0=(
                                float(y_fit_try[peak_idx]),
                                float(x_fit_try[peak_idx]),
                                float(bin_width_for_fit * 3),
                            ),
                            bounds=([0, -np.inf, 0], [np.inf, np.inf, np.inf]),
                        )

                        _params = fit_res_try["params"]
                        _params_err = fit_res_try["params_err"]
                        _has_nan = np.any(np.isnan(_params)) or np.any(np.isnan(_params_err))
                        if _has_nan:
                            _is_valid = False
                            log.debug(
                                "CI=%s: fit returned NaN params in bin %s",
                                ci, bin_label_dict[bin_idx]["label"],
                            )
                        else:
                            _sigma = abs(_params[2])
                            _sigma_rel_err_ok = _sigma <= 0 or _params_err[2] / _sigma < max_rel_err
                            _is_valid = _sigma_rel_err_ok
                            if not _sigma_rel_err_ok:
                                log.debug(
                                    "CI=%s: sigma rel. err. %.2f >= %.2f in bin %s",
                                    ci, _params_err[2] / _sigma, max_rel_err,
                                    bin_label_dict[bin_idx]["label"],
                                )
                        fit_candidates.append(
                            (ci, fit_res_try, x_fit_try, y_fit_try, y_fit_err_try, _is_valid)
                        )

                    # Select the candidate with the highest p-value (lowest chi2/ndf).
                    # Prefer valid candidates; fall back to best invalid if none are valid.
                    if fit_candidates:
                        valid_candidates = [c for c in fit_candidates if c[5]]
                        pool = valid_candidates if valid_candidates else fit_candidates
                        best = max(
                            pool,
                            key=lambda c: c[1]["p_value"] if not np.isnan(c[1]["p_value"]) else -np.inf,
                        )
                        ci_best, fit_res, x_fit, y_fit, y_fit_err, fit_accepted = best
                        log.debug(
                            "Selected CI=%s for bin %s (p_value=%.4f, chi2/ndf=%.2f/%d, accepted=%s)",
                            ci_best, bin_label_dict[bin_idx]["label"],
                            fit_res["p_value"], fit_res["chi2"], fit_res["dof"], fit_accepted,
                        )

                    if fit_res is not None:
                        # Exclude fit_func lambda — lambdas are not picklable and
                        # multiprocessing serializes the return value of this function.
                        gaussian_fits[bin_idx] = {
                            k: v for k, v in fit_res.items() if k != "fit_func"
                        }
                        if fit_candidates:
                            gaussian_fits[bin_idx]["ci_best"] = ci_best
                        if rebin_info_for_idx is not None:
                            gaussian_fits[bin_idx]["fit_rebin_info"] = rebin_info_for_idx
                        if fit_accepted:
                            resolutions[bin_idx] = fit_res["params"][2]
                            resolution_grid[bin_idx] = fit_res["params"][2]
                            resolutions_uncertainty[bin_idx] = fit_res["params_err"][2]
                            response_means[bin_idx] = fit_res["params"][1]
                            response_means_uncertainty[bin_idx] = fit_res["params_err"][1]
                        else:
                            log.warning(
                                "Fit quality too poor for bin %s for '%s' (best sigma rel. err. >= %.2f): setting resolution to NaN",
                                bin_label_dict[bin_idx]["label"], response_var_name, max_rel_err,
                            )
                            resolutions[bin_idx] = np.nan
                            resolution_grid[bin_idx] = np.nan
                            resolutions_uncertainty[bin_idx] = 0
                    else:
                        log.warning(
                            "Not enough valid points for Gaussian fit in bin %s for '%s': n_points=%d, total_counts=%.1f",
                            bin_label_dict[bin_idx]["label"],
                            response_var_name,
                            len(x_fit_base),
                            y_fit_base.sum(),
                        )
                        failed = True
                except Exception as e:
                    log.error("Failed to fit Gaussian for bin %s for '%s': %s", bin_label_dict[bin_idx]["label"], response_var_name, e)
                    failed = True
            else:
                try:
                    _ci_conf_level = cfg.get("ci_conf_level", 0.87)
                    resolution, lo_idx, hi_idx = Confidence_numpy(
                        hist=y_fit,
                        bins_mid=x_fit,
                        bin_width=bin_width_for_fit,
                        confLevel=_ci_conf_level,
                        return_bins=True,
                    )
                    resolutions[bin_idx] = resolution
                    resolution_grid[bin_idx] = resolution
                    ci_intervals[bin_idx] = (float(x_fit[lo_idx]), float(x_fit[hi_idx]))
                    ntot = np.sum(y_fit)
                    if ntot > 0:
                        # mean = sum(n_i * x_i) / N
                        mean = np.sum(y_fit * x_fit) / ntot
                        # mu2 = sum(n_i * (x_i - mean)^2) / N  [variance = RMS^2]
                        mu2 = np.sum(y_fit * (x_fit - mean) ** 2) / ntot
                        # mu4 = sum(n_i * (x_i - mean)^4) / N  [4th central moment]
                        mu4 = np.sum(y_fit * (x_fit - mean) ** 4) / ntot
                        # uncertainty on variance estimator (high N): sigma(s^2) = sqrt((mu4 - mu2^2) / N)
                        # uncertainty on RMS via error propagation d(RMS)/d(s^2) = 1/(2*RMS):
                        #   sigma(RMS) = sigma(s^2) / (2 * RMS)
                        #             = sqrt((mu4 - mu2^2) / N) / (2 * sqrt(mu2))
                        #             = sqrt((mu4 - mu2^2) / (4 * N * mu2))
                        # reduces to RMS/sqrt(2N) for a Gaussian where mu4 = 3*mu2^2
                        resolutions_uncertainty[bin_idx] = (
                            np.sqrt((mu4 - mu2**2) / (4 * ntot * mu2)) if mu2 > 0 else 0
                        )
                        response_means[bin_idx] = mean
                        # standard error of the mean = sigma / sqrt(N) = sqrt(mu2 / N)
                        response_means_uncertainty[bin_idx] = (
                            np.sqrt(mu2 / ntot) if mu2 > 0 else 0
                        )
                    else:
                        resolutions_uncertainty[bin_idx] = 0
                        response_means[bin_idx] = np.nan
                        response_means_uncertainty[bin_idx] = 0
                except Exception as e:
                    log.error("Failed to compute resolution for bin %s for '%s': %s", bin_label_dict[bin_idx]["label"], response_var_name, e)
                    failed = True
        else:
            log.debug(
                "Skipping bin %s: insufficient data (%d points)",
                bin_label_dict[bin_idx]["label"],
                int(np.sum(hist_counts)),
            )
            failed = True

        if failed:
            log.warning("Setting resolution to NaN for bin %s for '%s' due to failure", bin_label_dict[bin_idx]["label"], response_var_name)
            resolutions[bin_idx] = np.nan
            resolution_grid[bin_idx] = np.nan
            resolutions_uncertainty[bin_idx] = 0
            response_means.setdefault(bin_idx, np.nan)
            response_means_uncertainty.setdefault(bin_idx, 0)

    # Optionally override the per-bin response means with a configurable
    # computation method (Gaussian fit / histogram moments / Mean storage).
    # Must run before normalization so the normalization uses the chosen mean.
    _apply_response_mean_method(
        response_var_name=response_var_name,
        h_response=h_response,
        h_mean=h_mean,
        response_means=response_means,
        response_means_uncertainty=response_means_uncertainty,
        gaussian_fits=gaussian_fits,
        bin_shape=bin_shape,
        response_axis=response_axis,
        response_vars=response_vars,
    )

    # Optionally normalize the resolution by the mean of the response and
    # apply any of the supported per-bin uncertainty adjustments. See
    # _apply_response_normalization_and_extra_uncertainty for the full list.
    _apply_response_normalization_and_extra_uncertainty(
        response_var_name=response_var_name,
        resolutions=resolutions,
        resolutions_uncertainty=resolutions_uncertainty,
        resolution_grid=resolution_grid,
        response_means=response_means,
        response_means_uncertainty=response_means_uncertainty,
        response_vars=response_vars,
        gaussian_fits=gaussian_fits,
    )

    return (
        response_var_name,
        {
            "response_var": response_var_name,
            "bin_edges": bin_edges_dict,
            "bin_var_names": bin_var_names,
            "resolutions": resolutions,
            "resolution_grid": resolution_grid,
            "gaussian_fits": gaussian_fits,
            "rebin_infos": rebin_infos,
            "ci_intervals": ci_intervals,
            "resolutions_uncertainty": resolutions_uncertainty,
            "response_means": response_means,
            "response_means_uncertainty": response_means_uncertainty,
        },
    )


def compute_binned_resolution_from_histograms(
    h_dict, bin_var_names, response_vars, h_mean_dict=None
):
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
    h_mean_dict : dict, optional
        Optional mapping of response_var_name -> hist.Hist with Mean storage
        (axes = bin variables only). Used when a response variable's config
        sets ``response_mean_method="mean_storage"``.

    Returns
    -------
    dict
        Maps response_var_name -> resolution_result dict with keys:
            - "bin_edges": dict mapping bin_var_name to bin edges
            - "bin_var_names": list of bin variable names
            - "resolutions": dict of bin indices -> resolution values
            - "resolution_grid": ndarray with resolution values
    """
    h_mean_dict = h_mean_dict or {}
    # Prepare arguments for parallel processing
    args_list = [
        (
            response_var_name,
            h_response,
            bin_var_names,
            response_vars,
            h_mean_dict.get(response_var_name),
        )
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


def rebin_histogram(h, quantile_range=(0.01, 0.99), min_events_per_bin=20):
    """
    Rebin a 1D histogram automatically based on statistics.

    The target number of bins is derived from the total event count so that each
    bin has roughly ``min_events_per_bin`` entries, keeping statistical
    fluctuations small.

    Parameters
    ----------
    h : hist.Hist
        1D histogram to rebin.
    quantile_range : tuple of float
        Fraction of the distribution to keep (tails are trimmed).
    min_events_per_bin : int
        Desired minimum events per bin; drives the auto target_bins computation.

    Returns
    -------
    h_rebinned : hist.Hist
        Rebinned histogram.
    rebin_info : tuple or None
        ``(lo_edge, hi_edge, rebin_factor)`` describing what was applied, or
        ``None`` if the histogram was empty and no rebinning was done.
    """
    counts = h.values()
    total = counts.sum()

    if total == 0:
        log.warning("Histogram is empty, returning original without rebinning")
        return h, None

    target_bins = max(5, int(total / max(1, min_events_per_bin)))

    axis = h.axes[0]
    centers = axis.centers

    cdf = np.cumsum(counts) / total
    lo_idx = max(0, np.searchsorted(cdf, quantile_range[0]) - 1)
    hi_idx = min(len(centers), np.searchsorted(cdf, quantile_range[1]) + 1)
    n_bins_in_range = hi_idx - lo_idx

    # Floor division so the output has *at least* target_bins bins.
    rebin_factor = max(1, n_bins_in_range // target_bins)
    remainder = n_bins_in_range % rebin_factor
    if remainder != 0:
        hi_idx = min(len(centers), hi_idx + (rebin_factor - remainder))

    lo_edge = axis.edges[lo_idx]
    hi_edge = axis.edges[hi_idx]

    h_rebinned = h[lo_edge * 1j : hi_edge * 1j : rebin_factor * 1j]
    return h_rebinned, (lo_edge, hi_edge, rebin_factor)


def plot_variable_slices(
    h_dict,
    variables_dict,
    bin_var_configs,
    output_dir="",
    category="",
    year="",
    var_type="mapping",
    gaussian_fit_results=None,
    rebin_infos_results=None,
    resolution_results=None,
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

    # Group variables by their bin structure.
    # In mixed_mode, use name_plot of each bin axis as the key so that regular
    # and neutrino histograms (which share the same name_plot) are plotted
    # together. Canonical concrete axis names (from the first histogram in each
    # group) are stored separately for bin_var_configs lookups.
    variables_by_bin_name_plots = {}
    canonical_bin_axes = {}
    for var_name, h_var in h_dict.items():
        if var_name not in variables_dict:
            log.warning(
                "Variable '%s' found in histogram dict but not in variables config; skipping.",
                var_name,
            )
            continue
        n_bin_axes = h_var.ndim - 1
        concrete_axes = tuple(h_var.axes[i].name for i in range(n_bin_axes))
        if cfg["mixed_mode"]:
            group_key = tuple(
                bin_var_configs.get(ax, {}).get("name_plot", ax)
                for ax in concrete_axes
            )
        else:
            group_key = concrete_axes
        if group_key not in variables_by_bin_name_plots:
            variables_by_bin_name_plots[group_key] = {}
            canonical_bin_axes[group_key] = concrete_axes
        variables_by_bin_name_plots[group_key][var_name] = h_var

    plotters = []

    # Process each bin structure group separately
    for name_plot_key, group_dict in variables_by_bin_name_plots.items():
        bin_var_names = list(canonical_bin_axes[name_plot_key])
        n_bin_axes = len(bin_var_names)

        if n_bin_axes == 0:
            # No binning variables, just plot all variables in this group together
            series_dict = {}
            for var_name, h_var in group_dict.items():
                var_cfg = variables_dict[var_name]
                if np.sum(h_var.values()) == 0:
                    log.warning("Histogram %s is empty, skipping.", var_name)
                    continue
                series_dict[var_name] = {
                    "data": h_var,
                    "style": {
                        "color": get_color(var_name),
                        "legend_name": var_cfg.get("legend_name", var_name),
                    },
                }

            if not series_dict:
                log.warning(
                    "All histograms in group are empty for unbinned case, skipping plot."
                )
                continue

            # Create output filename
            output_name = f"{output_dir}/histo_{variables_dict[var_name]['name_plot']}_slice_{category}"

            log.debug(
                "Creating plot for unbinned variables: %s...", list(group_dict.keys())
            )
            # Plot using HEPPlotter
            plotter = (
                HEPPlotter()
                .set_plot_config(
                    cmstext="Simulation Preliminary" if not DP_NOTE_PLOTS else "Simulation\nPreliminary",
                    cmstext_font_size=35,
                    lumitext=f"{year} (13.6 TeV)" if year else "",
                    **({"cmstext_loc": 2} if DP_NOTE_PLOTS else {}),
                )
                .set_data(series_dict, plot_type="1d")
                .set_labels(
                    xlabel=f"{list(group_dict.values())[0].axes[-1].label}",
                    ylabel="Events",
                )
                .set_output(output_name)
                .set_options(grid=not DP_NOTE_PLOTS, legend=True)
            )
            plotters.append(plotter)
            continue

        # Get the first histogram in this group to determine bin shape
        first_h = list(group_dict.values())[0]
        bin_shape = tuple(len(first_h.axes[i]) for i in range(n_bin_axes))

        bin_edges_dict_plot = {
            first_h.axes[i].name: np.array(first_h.axes[i].edges)
            for i in range(n_bin_axes)
        }
        bin_label_dict = build_bin_label_dict(bin_shape, bin_var_names, bin_edges_dict_plot)
        # Iterate over all bin combinations
        for bin_idx in np.ndindex(bin_shape):
            # Build annotation text from bin ranges
            annotation_text = f"{PUPPI_JET_STRING}\n"
            for bin_var_name, (low_edge, high_edge) in bin_label_dict[bin_idx]["edges"].items():
                label = bin_var_configs[bin_var_name].get("label", bin_var_name)
                annotation_text += _format_bin_annotation_line(low_edge, high_edge, label) + "\n"

            # Collect histograms for this bin combination
            series_dict = {}
            for var_name, h_var in group_dict.items():
                var_cfg = variables_dict[var_name]

                if _bin_excluded_by_eta_max(
                    bin_idx, bin_var_names, bin_edges_dict_plot, var_cfg.get("eta_max")
                ):
                    continue

                # Extract 1D histogram for this bin combination
                slice_tuple = bin_idx + (slice(None),)
                h_1d = h_var[slice_tuple]

                if np.sum(h_1d.values()) == 0:
                    log.warning(
                        "Histogram %s is empty for bin combination %s, skipping.",
                        var_name,
                        bin_label_dict[bin_idx]["label"],
                    )
                    continue

                series_dict[var_name] = {
                    "data": h_1d,
                    "style": {
                        "color": get_color(var_name),
                        "legend_name": var_cfg.get("legend_name", var_name),
                        "linewidth": 2,
                    },
                }

            if not series_dict:
                log.warning(
                    "All histograms in group %s are empty for bin combination %s, skipping plot.",
                    bin_var_names,
                    bin_label_dict[bin_idx]["label"],
                )
                continue

            if variables_dict[var_name].get("rebin_for_plotting", False):
                rebinned_series_dict = {}
                for vn, var_data in series_dict.items():
                    var_fit = (
                        (gaussian_fit_results or {}).get(vn, {}).get(bin_idx)
                        if gaussian_fit_results is not None
                        else None
                    )
                    rebin_info = var_fit.get("fit_rebin_info") if var_fit is not None else None
                    if rebin_info is None and rebin_infos_results is not None:
                        rebin_info = (rebin_infos_results or {}).get(vn, {}).get(bin_idx)
                    if rebin_info is not None:
                        lo_edge, hi_edge, rf = rebin_info
                        h_rb = var_data["data"][lo_edge * 1j : hi_edge * 1j : rf * 1j]
                    else:
                        h_rb, _ = rebin_histogram(var_data["data"])
                    rebinned_series_dict[vn] = {"data": h_rb, "style": var_data["style"]}
                series_dict = rebinned_series_dict

            # Create output filename with bin ranges
            filename_parts = [
                f"histo_{variables_dict[var_name]['name_plot']}_slice",
                category,
            ]
            for bin_var_name, (low_edge, high_edge) in bin_label_dict[bin_idx]["edges"].items():
                low_edge_str = f"{low_edge}".replace(".", "p").replace("-", "m")
                high_edge_str = f"{high_edge}".replace(".", "p").replace("-", "m")
                filename_parts.append(
                    f"{bin_var_configs[bin_var_name]['name_plot']}_{low_edge_str}to{high_edge_str}"
                )

            output_name = f"{output_dir}/{'_'.join(filename_parts)}"

            # Get the axis label from tlineahe first variable in this group
            axis_label = list(group_dict.values())[0].axes[-1].label

            log.debug(
                "Creating plot for bin combination %s with variables: %s...",
                bin_label_dict[bin_idx]["label"],
                list(group_dict.keys()),
            )
            # Plot using HEPPlotter
            _ann_y_hist = 0.75 if DP_NOTE_PLOTS else 0.95
            plotter = (
                HEPPlotter()
                .set_plot_config(
                    cmstext="Simulation Preliminary" if not DP_NOTE_PLOTS else "Simulation\nPreliminary",
                    cmstext_font_size=30,
                    lumitext=f"{year} (13.6 TeV)" if year else "",
                    figsize=(13, 13) if not DP_NOTE_PLOTS else None,
                    **({"cmstext_loc": 2} if DP_NOTE_PLOTS else {}),
                )
                .set_data(series_dict, plot_type="1d")
                .set_labels(xlabel=axis_label, ylabel="Events")
                .set_output(output_name)
                .set_options(
                    grid=not DP_NOTE_PLOTS,
                    legend=True,
                    legend_loc="upper right",
                    y_log=not DP_NOTE_PLOTS,
                    split_legend=False,
                    ylim_bottom_value=1,
                    ylim_top_factor=1.8,
                )
                .add_annotation(
                    0.05,
                    _ann_y_hist,
                    annotation_text,
                    coord_type="axes",
                    verticalalignment="top",
                    fontsize=25,
                )
            )

            # Annotate sigma and overlay Gaussian fit curves when available
            if resolution_results is not None and not DP_NOTE_PLOTS:
                for i, var_name in enumerate(group_dict):
                    res_entry = (resolution_results.get(var_name) or {}).get(bin_idx)
                    if res_entry is None or np.isnan(res_entry[0]):
                        continue
                    sigma_val, sigma_err, ci_bounds = res_entry
                    col = i // 3
                    row = i % 3
                    ann_x = 0.05 + col * 0.3
                    ann_y = 0.75 - row * 0.07
                    sigma_label = f"$\\sigma={abs(sigma_val):.3f} \\pm {sigma_err:.3f}$"

                    if gaussian_fit_results is not None and cfg["gaussian_fit_resolution"]:
                        var_fits = gaussian_fit_results.get(var_name, {})
                        fit_res = var_fits.get(bin_idx)
                        if fit_res is not None and not np.any(np.isnan(fit_res["params"])):
                            sigma_label += (
                                f"\n$\\chi^2/\\mathrm{{ndf}}={fit_res['chi2']:.0f}/{fit_res['dof']}$, "
                                f"$p={fit_res['p_value']:.2f}$"
                            )
                            h_1d = series_dict[var_name]["data"]
                            axis = h_1d.axes[0]
                            x_fine = np.linspace(axis.edges[0], axis.edges[-1], 300)
                            y_fine = gaussian_model(
                                x_fine,
                                fit_res["params"][0],
                                fit_res["params"][1],
                                fit_res["params"][2],
                            )
                            plotter.add_curve(
                                x_fine,
                                y_fine,
                                color=get_color(var_name),
                                linestyle="--",
                                linewidth=2,
                            )
                            if cfg.get("gaussian_fit_cut_tails", False):
                                plotter.add_line(
                                    "v", x=fit_res["x_min"], color=get_color(var_name), linestyle="dotted"
                                )
                                plotter.add_line(
                                    "v", x=fit_res["x_max"], color=get_color(var_name), linestyle="dotted"
                                )
                                _ci_label = fit_res.get("ci_best")
                                if _ci_label is not None:
                                    _ci_str = f"CI={_ci_label:.2f}"
                                    _amp = fit_res["params"][0]
                                    for _x_line, _ha in ((fit_res["x_min"], "right"), (fit_res["x_max"], "left")):
                                        plotter.add_annotation(
                                            _x_line,
                                            _amp,
                                            _ci_str,
                                            coord_type="data",
                                            fontsize=10,
                                            color=get_color(var_name),
                                            ha=_ha,
                                            va="center",
                                            rotation=90,
                                        )

                    if ci_bounds is not None and not (
                        gaussian_fit_results is not None and cfg["gaussian_fit_resolution"]
                    ) and False:
                        x_lo, x_hi = ci_bounds
                        _ci_str = f"CI={cfg.get('ci_conf_level', 0.87):.2f}"
                        plotter.add_line(
                            "v", x=x_lo, color=get_color(var_name), linestyle="dotted"
                        )
                        plotter.add_line(
                            "v", x=x_hi, color=get_color(var_name), linestyle="dotted"
                        )
                        h_1d = series_dict[var_name]["data"]
                        _amp = float(h_1d.values().max()) if h_1d.values().max() > 0 else 1.0
                        for _x_line, _ha in ((x_lo, "right"), (x_hi, "left")):
                            plotter.add_annotation(
                                _x_line,
                                _amp,
                                _ci_str,
                                coord_type="data",
                                fontsize=10,
                                color=get_color(var_name),
                                ha=_ha,
                                va="center",
                                rotation=90,
                            )

                    plotter.add_annotation(
                        ann_x,
                        ann_y,
                        sigma_label,
                        coord_type="axes",
                        verticalalignment="top",
                        fontsize=15,
                        color=get_color(var_name),
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
        if mv not in mapping_vars:
            log.warning(
                "Mapping variable '%s' found in coffea file but not in mapping_vars config; skipping.",
                mv,
            )
            continue
        map_cfg = mapping_vars[mv]
        existing_axes = {ax.name for ax in h_mean.axes}
        active_bin_vars = [v for v in map_cfg["bin_vars"] if v in existing_axes]

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


def _make_mean_label(label):
    parts = label.split("$")
    math = parts[1] if len(parts) >= 3 else label.strip("$")
    units = parts[2].strip() if len(parts) >= 3 and parts[2].strip() else ""
    return f"$< {math} >$" + (f" {units}" if units else "")


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
    h_mean_bin_var_dict=None,
):
    log.info("Plotting linear fit for mapping variable '%s' vs its bin variable...", var_name)
    var_cfg = mapping_vars[var_name]
    x_var_name = var_cfg["bin_vars"][0]
    x_cfg = bin_var_configs[x_var_name]
    
    # Optionally use mean bin-var value per bin instead of the bin center.
    use_mean_bin_var = var_cfg.get("fit_with_mean_bin_var", False)

    x_label = x_cfg["label"]
    x_name = x_cfg["name_plot"]
    y_label = var_cfg["label"]
    y_name = var_cfg["name_plot"]

    y_mean_label = _make_mean_label(y_label)
    x_plot_label = _make_mean_label(x_label) if use_mean_bin_var else x_label

    # 2D histogram of mapping variable vs its bin variable
    var_2d = compute_projection(h_dict, mapping_vars, var_name)
    if var_2d.ndim < 2:
        log.warning("2D histogram for %s has only %d axis after squeezing; skipping 2D plot.", var_name, var_2d.ndim)
    elif np.sum(var_2d.values()) == 0:
        log.warning("2D histogram for %s is empty, skipping plot.", var_name)
    else:
        (
            HEPPlotter()
            .set_plot_config(cmstext="Simulation Preliminary",lumitext=f"{year} (13.6 TeV)")
            .set_options(legend=False, cbar_log=False)
            .set_output(f"{output_dir}/{y_name}_vs_{x_name}_2d_{category}")
            .set_data({"data": {"data": var_2d, "style": {}}}, plot_type="2d")
            .set_labels(x_label, y_label)
        ).run()

    x_values = None
    if use_mean_bin_var and h_mean_bin_var_dict and x_var_name in h_mean_bin_var_dict:
        bv_view = h_mean_bin_var_dict[x_var_name].view()
        x_values = np.where(bv_view.count > 0, bv_view.value, np.nan)
    elif use_mean_bin_var:
        log.warning(
            "fit_with_mean_bin_var requested for %s but mean bin-var histogram for '%s' "
            "is not available; falling back to bin centers.",
            var_name,
            x_var_name,
        )

    fit_results = perform_linear_fit(results[var_name]["mean"], x_values=x_values)
    log.debug("Linear fit results: %s", fit_results)

    x_axis = results[var_name]["mean"].axes[0]
    x_lin = np.linspace(x_values[0]*0.9, x_values[-1]*1.1, 100) if x_values is not None else np.linspace(x_axis.edges[0], x_axis.edges[-1], 100)
    x_centers = x_values if x_values is not None else (x_axis.edges[:-1] + x_axis.edges[1:]) / 2
    x_width = None if x_values is not None else (x_axis.widths / 2)
    mean_dict = {}
    mean_counts = results[var_name]["mean"].view().count
    if not np.all(mean_counts == 0):
        mean_dict["data"] = {
            "data": {
                "x": [x_centers, x_width],
                "y": [
                    results[var_name]["mean"].view().value,
                    np.sqrt(
                        results[var_name]["mean"].view().variance
                        / np.where(
                            mean_counts > 0,
                            mean_counts,
                            np.inf,
                        )
                    ),
                ],
            },
            "style": {"fmt": "o"},
        }
        mean_dict["linear fit"] = {
            "data": {
                "x": [x_lin, np.zeros(100)],
                "y": [fit_results["fit_func"](x_lin), np.zeros(100)],
            },
            "style": {"linestyle": "-", "fmt": ""},
        }

    if not mean_dict:
        log.warning(
            "Mean histogram for %s has no data, skipping linear fit plot.", var_name
        )
    else:
        (
            HEPPlotter()
            .set_plot_config(cmstext="Simulation Preliminary",lumitext=f"{year} (13.6 TeV)")
            .set_options(set_ylim=False)
            .set_output(f"{output_dir}/{y_name}_vs_{x_name}_fit_{category}")
            .set_data(mean_dict, plot_type="graph")
            .set_labels(x_plot_label, y_mean_label)
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
            nan_indices = [
                i
                for i, v in enumerate(row)
                if not isinstance(v, (int, np.integer)) and np.isnan(float(v))
            ]
            if nan_indices:
                log.error(
                    "NaN values detected in txt row at column(s) %s: %s — file '%s' may be corrupt",
                    nan_indices,
                    row,
                    output_path,
                )
            f.write("  ".join(_fmt(v) for v in row) + "\n")
    log.info("Saved resolution txt to %s", output_path)


def plot_1d_inclusive_distributions(
    *,
    response_h_dict,
    mapping_h_dict,
    available_response_vars,
    mapping_vars,
    bin_var_configs,
    output_dir,
    category,
    year,
):
    """
    Plot inclusive 1D distributions of response, mapping, and binning variables.

    For each variable, the ND histogram is projected onto the variable axis alone
    (all bin axes are summed over), giving a dataset-inclusive 1D distribution.
    Response and mapping variables that share the same ``name_plot`` are overlaid
    on a single canvas. Each binning variable gets its own plot.
    """
    os.makedirs(output_dir, exist_ok=True)
    log.info("Plotting inclusive 1D distributions for category '%s' in %s", category, output_dir)
    annotation_base = PUPPI_JET_STRING + "\nInclusive"

    def _make_plotter(series_dict, xlabel, output_name):
        _ann_y_incl = 0.75 if DP_NOTE_PLOTS else 0.95
        return (
            HEPPlotter()
            .set_plot_config(
                cmstext="Simulation Preliminary" if not DP_NOTE_PLOTS else "Simulation\nPreliminary",
                cmstext_font_size=30,
                lumitext=f"{year} (13.6 TeV)" if year else "",
                figsize=(13, 13),
                **({"cmstext_loc": 2} if DP_NOTE_PLOTS else {}),
            )
            .set_data(series_dict, plot_type="1d")
            .set_labels(xlabel=xlabel, ylabel="Events")
            .set_output(output_name)
            .set_options(
                grid=not DP_NOTE_PLOTS,
                legend=True,
                legend_loc="upper right",
                y_log=not DP_NOTE_PLOTS,
                split_legend=False,
                ylim_bottom_value=1,
            )
            .add_annotation(
                0.05,
                _ann_y_incl,
                annotation_base,
                coord_type="axes",
                verticalalignment="top",
                fontsize=20,
            )
        )

    # --- response variables: group by name_plot, overlay on one canvas ---
    response_groups = {}
    for var_name, h_nd in response_h_dict.items():
        if var_name not in available_response_vars:
            continue
        name_plot = available_response_vars[var_name].get("name_plot", var_name)
        response_groups.setdefault(name_plot, {})[var_name] = h_nd

    for name_plot, group in response_groups.items():
        series_dict = {}
        xlabel = None
        for var_name, h_nd in group.items():
            h_1d = h_nd.project(var_name)
            if np.sum(h_1d.values()) == 0:
                continue
            if available_response_vars[var_name].get("rebin_for_plotting", False):
                h_1d, _ = rebin_histogram(h_1d)
            if xlabel is None:
                xlabel = h_nd.axes[var_name].label
            series_dict[var_name] = {
                "data": h_1d,
                "style": {
                    "color": get_color(var_name),
                    "legend_name": available_response_vars[var_name].get("legend_name", var_name),
                    "linewidth": 2,
                },
            }
        if not series_dict:
            continue
        output_name = f"{output_dir}/inclusive_response_{name_plot}_{category}"
        _make_plotter(series_dict, xlabel or name_plot, output_name).run()

    # --- mapping variables: group by name_plot, overlay on one canvas ---
    mapping_groups = {}
    for var_name, h_nd in mapping_h_dict.items():
        if var_name not in mapping_vars:
            continue
        name_plot = mapping_vars[var_name].get("name_plot", var_name)
        mapping_groups.setdefault(name_plot, {})[var_name] = h_nd

    for name_plot, group in mapping_groups.items():
        series_dict = {}
        xlabel = None
        for var_name, h_nd in group.items():
            h_1d = h_nd.project(var_name)
            if np.sum(h_1d.values()) == 0:
                continue
            if mapping_vars[var_name].get("rebin_for_plotting", False):
                h_1d, _ = rebin_histogram(h_1d)
            if xlabel is None:
                xlabel = h_nd.axes[var_name].label
            series_dict[var_name] = {
                "data": h_1d,
                "style": {
                    "color": get_color(var_name),
                    "legend_name": mapping_vars[var_name].get("legend_name", var_name),
                    "linewidth": 2,
                },
            }
        if not series_dict:
            continue
        output_name = f"{output_dir}/inclusive_mapping_{name_plot}_{category}"
        _make_plotter(series_dict, xlabel or name_plot, output_name).run()

    # --- bin variables: project from any available ND histogram ---
    seen_bin_vars = {}
    for h_nd in list(mapping_h_dict.values()) + list(response_h_dict.values()):
        for axis in h_nd.axes:
            bv_name = axis.name
            if bv_name in bin_var_configs and bv_name not in seen_bin_vars:
                h_1d = h_nd.project(bv_name)
                if np.sum(h_1d.values()) > 0:
                    seen_bin_vars[bv_name] = h_1d

    for bv_name, h_1d in seen_bin_vars.items():
        bv_cfg = bin_var_configs[bv_name]
        xlabel = bv_cfg.get("label", bv_name)
        series_dict = {
            bv_name: {
                "data": h_1d,
                "style": {
                    "color": None,
                    "legend_name": xlabel,
                    "linewidth": 2,
                },
            }
        }
        output_name = f"{output_dir}/inclusive_binvar_{bv_cfg['name_plot']}_{category}"
        _make_plotter(series_dict, xlabel, output_name).run()


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

        rebin_infos_results = (
            {
                var_name: result.get("rebin_infos", {})
                for var_name, result in response_tot_dict.items()
            }
            if response_tot_dict
            else None
        )

        resolution_results = (
            {
                var_name: {
                    bin_idx: (
                        res,
                        result.get("resolutions_uncertainty", {}).get(bin_idx, 0),
                        result.get("ci_intervals", {}).get(bin_idx),
                    )
                    for bin_idx, res in result.get("resolutions", {}).items()
                }
                for var_name, result in response_tot_dict.items()
            }
            if response_tot_dict
            else None
        )

        plot_variable_slices(
            h_dict=response_h_dict,
            variables_dict=available_response_vars,
            bin_var_configs=bin_vars,
            output_dir=response_histogram_dir,
            category=category,
            year=year,
            var_type="response",
            gaussian_fit_results=gaussian_fit_results,
            rebin_infos_results=rebin_infos_results,
            resolution_results=resolution_results,
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
        h_mean_bin_var_dict = cat_data.get("h_mean_bin_var_dict", {})

        # Load linear fit results from file if available, otherwise recompute.
        # Skip entirely when --refit is active: the refit block below will
        # redo the fits after merging histogram bins.
        if "linear_fit_maps" in cat_data and not args.refit:
            linear_fit_maps = cat_data["linear_fit_maps"]
            mapped_bin_edges = cat_data.get("mapped_bin_edges", {})
        elif not args.refit:
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
                    h_mean_bin_var_dict=h_mean_bin_var_dict,
                )
                mapped_bin_edges[lf_var_name] = linear_fit_maps[lf_var_name][
                    "fit_func"
                ](cfg["bin_variables"][lf_var_cfg["bin_vars"][0]]["bin_edges"])

        # When --refit is requested, merge histogram bins from the YAML config and
        # recompute results / profile means so that the new binning is fully consistent.
        # NOTE: plot_mapping_variable_histograms is called *after* this block so that it
        # always receives an h_dict whose axes match cfg["bin_variables_mixed"].
        if args.refit:
            log.info("--refit: merging histogram bins from YAML config and recomputing profile means...")
            bin_var_configs_all = cfg["bin_variables_mixed"]

            # Merge mapping-variable Count histograms so linear fits stay consistent.
            # Single-bin axes (after merging) are squeezed out so they don't appear
            # in plot filenames or annotations.
            h_dict = {
                vn: squeeze_single_bin_axes(merge_nd_histogram_bins(h, bin_var_configs_all)[0], bin_var_configs_all)[0]
                for vn, h in h_dict.items()
            }

            # Recompute profile means from the merged Mean histograms (if saved).
            if "h_mean_dict" in cat_data:
                h_mean_dict_merged = {
                    vn: squeeze_single_bin_axes(merge_nd_histogram_bins(h, bin_var_configs_all)[0], bin_var_configs_all)[0]
                    for vn, h in cat_data["h_mean_dict"].items()
                }
                results = profile_means(h_mean_dict_merged, cfg["mapping_variables"])
            else:
                log.warning(
                    "--refit: h_mean_dict not found in coffea file; "
                    "profile means for the original binning will be reused."
                )

            # Merge mean bin-var histograms if present.
            if h_mean_bin_var_dict:
                h_mean_bin_var_dict = {
                    vn: squeeze_single_bin_axes(merge_nd_histogram_bins(h, bin_var_configs_all)[0], bin_var_configs_all)[0]
                    for vn, h in h_mean_bin_var_dict.items()
                }

            # Redo linear fits with the merged histograms.
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
                    h_mean_bin_var_dict=h_mean_bin_var_dict,
                )
                mapped_bin_edges[lf_var_name] = linear_fit_maps[lf_var_name]["fit_func"](
                    cfg["bin_variables"][lf_var_cfg["bin_vars"][0]]["bin_edges"]
                )

        plot_mapping_variable_histograms(category=category, year=year, h_dict=h_dict)

        # Process response variables from loaded data
        for mode in ["regular", "neutrino", "mixed"]:
            if mode not in cat_data["response_data"]:
                continue

            resp_data = cat_data["response_data"][mode]
            response_h_dict = resp_data["response_h_dict"]
            response_h_mean_dict = resp_data.get("response_h_mean_dict", {}) or {}
            available_response_vars = resp_data["available_response_vars"]

            bin_vars = get_bin_vars_for_mode(mode)

            if args.refit:
                log.info("--refit mode='%s': merging response histogram bins...", mode)
                response_h_dict = {
                    vn: squeeze_single_bin_axes(merge_nd_histogram_bins(h, bin_vars)[0], bin_vars)[0]
                    for vn, h in response_h_dict.items()
                }
                if response_h_mean_dict:
                    response_h_mean_dict = {
                        vn: squeeze_single_bin_axes(merge_nd_histogram_bins(h, bin_vars)[0], bin_vars)[0]
                        for vn, h in response_h_mean_dict.items()
                    }
                log.info("--refit mode='%s': recomputing Gaussian fit...", mode)
                response_tot_dict = compute_binned_resolution_from_histograms(
                    h_dict=response_h_dict,
                    bin_var_names=list(bin_vars.keys()),
                    response_vars=available_response_vars,
                    h_mean_dict=response_h_mean_dict,
                )
            else:
                response_tot_dict = resp_data["response_tot_dict"]

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
            inclusive_dir = f"{args.output}/inclusive_1d_{category}_{mode}"
            plot_1d_inclusive_distributions(
                response_h_dict=response_h_dict,
                mapping_h_dict=h_dict,
                available_response_vars=available_response_vars,
                mapping_vars=cfg["mapping_variables"],
                bin_var_configs=bin_vars,
                output_dir=inclusive_dir,
                category=category,
                year=year,
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
        h_dict, h_mean_dict, h_mean_bin_var_dict = create_ND_histo(
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
                h_mean_bin_var_dict=h_mean_bin_var_dict,
            )
            mapped_bin_edges[lf_var_name] = linear_fit_maps[lf_var_name]["fit_func"](
                cfg["bin_variables"][lf_var_cfg["bin_vars"][0]]["bin_edges"]
            )

        # Store category data for saving
        category_data = {
            "h_dict": h_dict,
            "h_mean_dict": h_mean_dict,
            "h_mean_bin_var_dict": h_mean_bin_var_dict,
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
            response_h_dict, response_h_mean_dict, _ = create_ND_histo(
                variables_dict=available_response_vars,
                data=col_var_flatten,
                bin_var_configs=bin_vars,
            )

            # Compute resolutions from the ND histograms
            response_tot_dict = compute_binned_resolution_from_histograms(
                h_dict=response_h_dict,
                bin_var_names=bin_var_names,
                response_vars=available_response_vars,
                h_mean_dict=response_h_mean_dict,
            )

            # Store response data for this category
            response_type_key = mode
            category_data["response_data"][response_type_key] = {
                "response_h_dict": response_h_dict,
                "response_h_mean_dict": response_h_mean_dict,
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
                mode=response_type,
                response_h_dict=resp_data["response_h_dict"],
                response_tot_dict=resp_data["response_tot_dict"],
                available_response_vars=resp_data["available_response_vars"],
                bin_vars=resp_data["bin_vars"],
                results_for_mapping=results,
                linear_fit_maps=linear_fit_maps,
            )
            inclusive_dir = f"{args.output}/inclusive_1d_{category}_{response_type}"
            plot_1d_inclusive_distributions(
                response_h_dict=resp_data["response_h_dict"],
                mapping_h_dict=h_dict,
                available_response_vars=resp_data["available_response_vars"],
                mapping_vars=cfg["mapping_variables"],
                bin_var_configs=resp_data["bin_vars"],
                output_dir=inclusive_dir,
                category=category,
                year=year,
            )

    log.info("All done!")


if __name__ == "__main__":

    main()
