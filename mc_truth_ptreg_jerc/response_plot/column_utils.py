"""Shared helpers to turn dumped columns (parquet) into binned results.

Both the JER MC study (``plot_jer_mc.py``) and the MC truth derivation
(``response_mc_truth.py``) follow the same pattern:

1. read the flat columns produced by PocketCoffea (``dump_columns_as_arrays_per_chunk``),
2. flatten them to one entry per jet,
3. fill N-dimensional histograms whose last axis is the quantity of interest
   (the response) and whose other axes are the binning variables,
4. extract one number per bin (a median for the MC truth, a width for the JER),
5. fit and plot.

This module collects the generic pieces of steps 1-4 so that a new analysis
only has to describe its binning in YAML. ``response_mc_truth.py`` is built
entirely on top of it.
"""

import copy
import logging
import os

import hist
import numpy as np
import yaml
from scipy.optimize import curve_fit
from scipy.stats import chi2 as chi2_dist

from mc_truth_ptreg_jerc.response_plot.confidence import Confidence_numpy
from mc_truth_ptreg_jerc.response_plot.pol_functions import pol_functions_dict

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# logging and configuration
# ---------------------------------------------------------------------------
def setup_logging(output_dir, log_name, loggers):
    """Attach a shared console (INFO) and file (DEBUG) handler to *loggers*."""
    if not isinstance(loggers, (list, tuple)):
        loggers = [loggers]

    fmt = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    console.setFormatter(fmt)

    log_path = os.path.join(output_dir, log_name)
    file_handler = logging.FileHandler(log_path, mode="a")
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(fmt)

    for logger in loggers:
        logger.setLevel(logging.DEBUG)
        logger.handlers = [console, file_handler]
        logger.propagate = False

    loggers[0].info("Logging to %s", log_path)
    return log_path


def deep_merge(base, override):
    """Recursively merge *override* into *base*; lists are replaced, not extended."""
    result = dict(base)
    for key, value in override.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            result[key] = deep_merge(result[key], value)
        else:
            result[key] = value
    return result


def load_yaml_config(config_path=None, default_path=None, test_override=False):
    """Load the default YAML config and merge a user config on top of it.

    Scalars and lists from *config_path* replace the defaults; nested dicts
    (``bin_variables``, ``response_variables``, ...) are merged key by key, so a
    custom config can tweak a single option of a single variable.
    """
    if default_path is not None:
        with open(default_path) as f:
            cfg = yaml.safe_load(f)
        if config_path is not None:
            with open(config_path) as f:
                override = yaml.safe_load(f) or {}
            cfg = deep_merge(cfg, override)
    else:
        with open(config_path) as f:
            cfg = yaml.safe_load(f)

    if test_override:
        cfg["test"] = True

    if cfg.get("test", False):
        for key, test_key in (
            ("jet_eta_bins", "test_jet_eta_bins"),
            ("jet_pt_bins", "test_jet_pt_bins"),
        ):
            if test_key in cfg:
                cfg[key] = cfg[test_key]
        for var_cfg in cfg.get("bin_variables", {}).values():
            name = var_cfg.get("name_plot", "")
            if "eta" in name and "test_jet_eta_bins" in cfg:
                var_cfg["bin_edges"] = cfg["test_jet_eta_bins"]
            elif "pt" in name and "test_jet_pt_bins" in cfg:
                var_cfg["bin_edges"] = cfg["test_jet_pt_bins"]

    for var_cfg in cfg.get("bin_variables", {}).values():
        if "bin_edges" in var_cfg:
            var_cfg["bin_edges"] = np.array(var_cfg["bin_edges"], dtype=float)

    return cfg


# ---------------------------------------------------------------------------
# columns
# ---------------------------------------------------------------------------
def flatten_data(data):
    """Flatten every column to one entry per jet.

    Jet level columns (arrays of arrays) are concatenated; event level columns
    (flat arrays) are repeated once per jet in the event.

    Note
    ----
    When the input holds several jet collections with different multiplicities
    (``MatchedJets`` and ``MatchedJetsNeutrino``), each of them keeps its own
    length, but the event level columns are broadcast with the multiplicity of
    the *first* jet collection found. They can therefore only be used as binning
    variables together with that collection. The MC truth derivation only bins
    in jet level quantities (flavour, eta, pT), so this never bites there.
    """

    def is_ragged(arr):
        return arr.dtype == object and len(arr) > 0 and isinstance(arr[0], np.ndarray)

    jet_level_key = next((k for k, a in data.items() if is_ragged(a)), None)
    if jet_level_key is None:
        raise ValueError("No jet level (jagged) column found: nothing to flatten.")

    n_jets = np.array([len(jets) for jets in data[jet_level_key]])

    flat_data = {}
    for key, arr in data.items():
        if is_ragged(arr):
            flat_data[key] = np.concatenate(arr) if len(arr) else np.array([])
        else:
            flat_data[key] = np.repeat(arr, n_jets)
    return flat_data


def flavour_group_names(flavour_cfg):
    """Return the ordered list of flavour groups, including the catch-all."""
    names = list(flavour_cfg["groups"].keys())
    names.append(flavour_cfg.get("other_name", "other"))
    return names


def add_flavour_index(data, flavour_cfg, collections):
    """Add an integer flavour index column for every jet collection.

    The index points into :func:`flavour_group_names`; jets whose flavour is not
    listed in any group land in the catch-all group, so that summing the flavour
    axis of a histogram always gives back the inclusive sample.

    Returns the list of the column names that were added.
    """
    groups = flavour_cfg["groups"]
    field = flavour_cfg.get("field", "partonFlavour")
    use_abs = flavour_cfg.get("use_abs", True)
    other_index = len(groups)

    added = []
    for collection in collections:
        source = f"{collection}_{field}"
        target = f"{collection}_FlavourIndex"
        if source not in data:
            log.warning(
                "Flavour column '%s' not found: no flavour splitting for '%s'.",
                source,
                collection,
            )
            continue
        values = np.asarray(data[source])
        values = np.abs(values) if use_abs else values
        index = np.full(len(values), other_index, dtype=np.int64)
        for i, pdg_ids in enumerate(groups.values()):
            pdg_ids = [pdg_ids] if np.isscalar(pdg_ids) else list(pdg_ids)
            index[np.isin(values, pdg_ids)] = i
        data[target] = index
        added.append(target)
    return added


# ---------------------------------------------------------------------------
# histograms
# ---------------------------------------------------------------------------
def make_bin_axis(bin_var_name, bin_var_cfg):
    """Build a hist axis from a bin-variable configuration entry.

    ``categories`` in the configuration selects an integer categorical axis
    (used for the flavour splitting); otherwise a variable-width axis is built
    from ``bin_edges``.
    """
    if "categories" in bin_var_cfg:
        return hist.axis.IntCategory(
            list(range(len(bin_var_cfg["categories"]))),
            name=bin_var_name,
            label=bin_var_cfg.get("label", bin_var_name),
        )
    return hist.axis.Variable(
        bin_var_cfg["bin_edges"],
        name=bin_var_name,
        label=bin_var_cfg.get("label", bin_var_name),
        flow=False,
    )


def create_ND_histo(variables_dict, data, bin_var_configs):
    """Fill one N-D histogram (and one Mean histogram) per variable.

    Returns ``(h_dict, h_mean_dict)`` where ``h_dict[var]`` has axes
    ``[*bin_vars, var]`` with Count storage and ``h_mean_dict[var]`` has axes
    ``[*bin_vars]`` with Mean storage (the un-binned mean of ``var`` per bin).
    """
    h_dict = {}
    h_mean_dict = {}

    for i, (var_name, var_cfg) in enumerate(variables_dict.items(), start=1):
        if var_name not in data:
            log.warning("Variable '%s' not found in the columns, skipping.", var_name)
            continue

        bin_var_names = var_cfg.get("bin_vars", [])
        axes_bin = []
        fill_kwargs_bin = {}
        missing = False
        for bin_var_name in bin_var_names:
            if bin_var_name not in bin_var_configs:
                raise ValueError(
                    f"Bin variable '{bin_var_name}' is not defined in bin_variables"
                )
            if bin_var_name not in data:
                log.warning(
                    "Bin variable '%s' not found in the columns: skipping '%s'.",
                    bin_var_name,
                    var_name,
                )
                missing = True
                break
            axes_bin.append(make_bin_axis(bin_var_name, bin_var_configs[bin_var_name]))
            fill_kwargs_bin[bin_var_name] = data[bin_var_name]
        if missing:
            continue

        values = np.asarray(data[var_name], dtype=float)

        if var_cfg.get("mean_only", False):
            # only the per-bin mean is needed: skip the (potentially huge)
            # count histogram along the variable axis
            if not axes_bin:
                raise ValueError(
                    f"'{var_name}' is declared mean_only but has no binning axes"
                )
            log.info("Filling mean histogram %d/%d: '%s'", i, len(variables_dict), var_name)
            h_mean = hist.Hist(*axes_bin, storage=hist.storage.Mean())
            h_mean.fill(**fill_kwargs_bin, sample=values)
            h_mean_dict[var_name] = h_mean
            continue

        if "bin_limits" in var_cfg:
            low, high = var_cfg["bin_limits"]
        else:
            low, high = np.nanmin(values), np.nanmax(values)
            padding = (high - low) * 0.001
            low, high = low - padding, high + padding

        var_axis = hist.axis.Regular(
            var_cfg.get("N_bins", 100),
            low,
            high,
            name=var_name,
            label=var_cfg.get("label", var_name),
            flow=False,
        )

        log.info("Filling histogram %d/%d: '%s'", i, len(variables_dict), var_name)
        h = hist.Hist(*axes_bin, var_axis)
        fill_kwargs = copy.deepcopy(fill_kwargs_bin)
        fill_kwargs[var_name] = values
        h.fill(**fill_kwargs)
        h_dict[var_name] = h

        if axes_bin:
            h_mean = hist.Hist(*axes_bin, storage=hist.storage.Mean())
            h_mean.fill(**fill_kwargs_bin, sample=values)
            h_mean_dict[var_name] = h_mean

    return h_dict, h_mean_dict


def rebin_histogram(h, quantile_range=(0.01, 0.99), min_events_per_bin=20,
                    max_rebin_factor=None):
    """Rebin a 1D histogram so that each bin holds ~*min_events_per_bin* entries."""
    counts = h.values()
    total = counts.sum()
    if total == 0:
        return h, None

    target_bins = max(5, int(total / max(1, min_events_per_bin)))
    axis = h.axes[0]
    centers = axis.centers

    cdf = np.cumsum(counts) / total
    lo_idx = max(0, np.searchsorted(cdf, quantile_range[0]) - 1)
    hi_idx = min(len(centers), np.searchsorted(cdf, quantile_range[1]) + 1)
    n_bins_in_range = hi_idx - lo_idx

    rebin_factor = max(1, n_bins_in_range // target_bins)
    if max_rebin_factor is not None:
        rebin_factor = min(rebin_factor, max_rebin_factor)
    remainder = n_bins_in_range % rebin_factor
    if remainder != 0:
        hi_idx = min(len(centers), hi_idx + (rebin_factor - remainder))

    lo_edge = axis.edges[lo_idx]
    hi_edge = axis.edges[hi_idx]
    return h[lo_edge * 1j : hi_edge * 1j : rebin_factor * 1j], (
        lo_edge,
        hi_edge,
        rebin_factor,
    )


def build_bin_label_dict(bin_shape, bin_var_names, bin_edges_dict):
    """Map every bin index tuple to a readable label and its per-variable edges."""
    result = {}
    for bin_idx in np.ndindex(bin_shape):
        parts = []
        edges_per_var = {}
        for i, var_name in enumerate(bin_var_names):
            edges = bin_edges_dict[var_name]
            low, high = edges[bin_idx[i]], edges[bin_idx[i] + 1]
            edges_per_var[var_name] = (low, high)
            parts.append(f"{var_name}: [{low:g}, {high:g})")
        result[bin_idx] = {"label": ", ".join(parts), "edges": edges_per_var}
    return result


def format_bin_annotation_line(low_edge, high_edge, label):
    """Format one bin range as a single annotation line."""
    if float(low_edge).is_integer() and float(high_edge).is_integer():
        return f"{int(low_edge)} < {label} < {int(high_edge)}"
    return f"{low_edge:.3g} < {label} < {high_edge:.3g}"


def edge_to_string(edge):
    """Turn a bin edge into a file-name friendly string (-1.3 -> m1p3)."""
    if float(edge).is_integer():
        text = f"{int(edge)}"
    else:
        text = f"{edge}".replace(".", "p")
    return text.replace("-", "m")


# ---------------------------------------------------------------------------
# per-bin estimators
# ---------------------------------------------------------------------------
def quantile_from_counts(counts, edges, quantile):
    """Quantile of a binned distribution, interpolating inside the bin.

    Taking the centre of the bin in which the cumulative distribution crosses
    *quantile* quantizes the result to the bin width. Interpolating the
    cumulative distribution linearly across that bin removes the quantization,
    which lets the response axis stay coarse enough to fit in memory while
    keeping the median accurate to well below the bin width.
    """
    total = counts.sum()
    if total <= 0:
        return np.nan

    cdf = np.cumsum(counts) / total
    index = int(np.searchsorted(cdf, quantile, side="left"))
    index = min(index, len(counts) - 1)

    cdf_low = cdf[index - 1] if index > 0 else 0.0
    fraction_in_bin = counts[index] / total
    if fraction_in_bin <= 0:
        return edges[index]
    return edges[index] + (quantile - cdf_low) / fraction_in_bin * (
        edges[index + 1] - edges[index]
    )


def median_and_error(counts, edges):
    """Median of a binned distribution and its uncertainty.

    The uncertainty follows the usual estimate for the median of a
    distribution, ``1.253 * RMS / sqrt(N)``.
    """
    counts = np.asarray(counts, dtype=float)
    total = counts.sum()
    if total <= 0:
        return np.nan, np.nan

    median = quantile_from_counts(counts, edges, 0.5)

    centers = (np.asarray(edges[:-1]) + np.asarray(edges[1:])) / 2
    mean = np.average(centers, weights=counts)
    rms = np.sqrt(np.average((centers - mean) ** 2, weights=counts))
    return median, 1.253 * rms / np.sqrt(total)


def quantile_width(counts, edges, q_low=0.16, q_high=0.84):
    """Half the distance between two quantiles, i.e. the 68% half-width."""
    counts = np.asarray(counts, dtype=float)
    total = counts.sum()
    if total <= 0:
        return np.nan, np.nan

    low = quantile_from_counts(counts, edges, q_low)
    high = quantile_from_counts(counts, edges, q_high)
    width = (high - low) / 2
    # For a Gaussian the uncertainty on the width is sigma / sqrt(2N).
    return width, width / np.sqrt(2 * total)


def confidence_width(counts, centers, conf_level=0.87):
    """Width of the smallest interval containing *conf_level* of the entries."""
    total = counts.sum()
    if total <= 0 or len(centers) < 2:
        return np.nan, np.nan
    bin_width = centers[1] - centers[0]
    width = Confidence_numpy(
        np.asarray(counts, dtype=float), np.asarray(centers, dtype=float), bin_width,
        confLevel=conf_level,
    )
    return width, width / np.sqrt(2 * total)


def mean_from_mean_storage(h_mean, bin_idx):
    """Read value and uncertainty of the mean from a Mean-storage histogram."""
    view = h_mean.view()
    value = float(np.asarray(view.value)[bin_idx])
    count = float(np.asarray(view.count)[bin_idx])
    variance = float(np.asarray(view.variance)[bin_idx])
    if count <= 0:
        return np.nan, np.nan
    return value, np.sqrt(max(variance, 0) / count)


# ---------------------------------------------------------------------------
# fits
# ---------------------------------------------------------------------------
def perform_fit(x, y, y_err, fit_function, n_params, p0=None, maxfev=100000):
    """Weighted least-squares fit returning parameters and fit quality."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    invalid = {
        "params": np.full(n_params, np.nan),
        "params_err": np.full(n_params, np.nan),
        "chi2": np.nan,
        "ndof": 0,
        "p_value": np.nan,
        "x_min": np.nan,
        "x_max": np.nan,
    }
    if len(x) < n_params:
        return invalid

    use_sigma = y_err is not None and np.all(np.asarray(y_err) > 0)
    kwargs = {"p0": p0, "maxfev": maxfev}
    if use_sigma:
        kwargs["sigma"] = np.asarray(y_err, dtype=float)
        kwargs["absolute_sigma"] = True

    try:
        popt, pcov = curve_fit(fit_function, x, y, **kwargs)
    except (RuntimeError, ValueError) as exc:
        log.warning("Fit failed: %s", exc)
        return invalid

    params = np.asarray(popt, dtype=float)
    params_err = (
        np.sqrt(np.diag(pcov))
        if pcov is not None and np.size(pcov)
        else np.full_like(params, np.nan)
    )

    residuals = y - fit_function(x, *params)
    chi2 = (
        float(np.sum((residuals / np.asarray(y_err, dtype=float)) ** 2))
        if use_sigma
        else float(np.sum(residuals**2))
    )
    ndof = max(len(x) - n_params, 1)
    return {
        "params": params,
        "params_err": params_err,
        "chi2": chi2,
        "ndof": ndof,
        "p_value": float(chi2_dist.sf(chi2, ndof)),
        "x_min": float(np.min(x)),
        "x_max": float(np.max(x)),
    }


def fit_pol_log10(x, y, y_err, max_order=10, p_value_threshold=0.05):
    """Fit ``y`` versus ``log10(x)`` with polynomials of increasing order.

    The first order whose p-value exceeds *p_value_threshold* is accepted. If no
    order converges to an acceptable p-value the one with the largest p-value is
    kept (or, when all p-values are numerically zero, the one with the smallest
    chi-square). This is the same strategy as the original ``response.py``.

    Returns ``None`` when no fit could be performed at all.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    y_err = np.asarray(y_err, dtype=float)

    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(y_err) & (x > 0)
    x, y, y_err = x[mask], y[mask], y_err[mask]
    if len(x) < 2:
        return None

    candidates = []
    for order, pol_function in pol_functions_dict.items():
        if order > max_order or order + 1 >= len(x):
            break
        fit_result = perform_fit(
            x, y, y_err, pol_function, n_params=order + 1, p0=[1.0] * (order + 1)
        )
        if not np.all(np.isfinite(fit_result["params"])):
            continue
        fit_result["pol"] = order
        candidates.append(fit_result)
        if fit_result["p_value"] > p_value_threshold:
            break

    if not candidates:
        return None

    p_values = [c["p_value"] for c in candidates]
    best_p_value = np.nanmax(p_values) if np.any(np.isfinite(p_values)) else np.nan
    if np.isfinite(best_p_value) and best_p_value > 1e-7:
        best = candidates[int(np.nanargmax(p_values))]
    else:
        best = candidates[int(np.nanargmin([c["chi2"] for c in candidates]))]

    best = dict(best)
    best["x"] = x.tolist()
    best["y"] = y.tolist()
    best["y_err"] = y_err.tolist()
    best["jet_pt"] = [float(x[0]), float(x[-1])]
    best["parameters"] = np.asarray(best["params"]).tolist()
    best["errors"] = np.asarray(best["params_err"]).tolist()
    return best


def evaluate_pol(parameters, x):
    """Evaluate the polynomial in ``log10(x)`` described by *parameters*."""
    order = len(parameters) - 1
    return pol_functions_dict[order](np.asarray(x, dtype=float), *parameters)
