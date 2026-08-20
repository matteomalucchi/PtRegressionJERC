#!/usr/bin/env python3
"""Derive the MC truth (L2Relative) corrections from the dumped columns.

This is the column-based replacement of ``response.py``. It reads the parquet
files produced by ``mc_truth_columns_config.py`` and, in a single pass:

1. fills one N-dimensional histogram per response variable, with axes
   ``(flavour, eta^reco, pT^ptcl, response)`` -- so ParticleNet, UParT and their
   neutrino-inclusive versions are all treated together and no job has to be
   split by eta region or by flavour any more;
2. extracts, per bin, the median of the response, its inverse (the correction),
   the resolution and the mean reconstructed pT;
3. fits the inverse median as a polynomial in ``log10(pT^reco)`` and writes the
   JME ``L2Relative`` txt files, inclusively and per flavour;
4. produces the summary plots (median, inverse median, resolution) overlaying
   the corrections for a given flavour and overlaying the flavours for a given
   correction, plus the response histograms per bin.

It shares its structure -- and, through ``column_utils.py``, most of its
machinery -- with ``plot_jer_mc.py``: same YAML-driven binning, same
``HEPPlotter`` output, same ``--load`` mechanism to replot without re-reading
the parquet files.

Examples
--------
    python response_mc_truth.py -i /work/$USER/out_mc_truth_2023_preBPix \
        -o /work/$USER/out_mc_truth_2023_preBPix/mc_truth_plots -w 16

    python response_mc_truth.py -l <output>/mc_truth_results_baseline.coffea \
        -o <output>_replot
"""

import logging

logging.getLogger("matplotlib").setLevel(logging.WARNING)
logging.getLogger("boost_histogram").setLevel(logging.WARNING)

log = logging.getLogger("response_mc_truth")

import argparse
import json
import os
import shutil
from multiprocessing import Pool

import numpy as np
from coffea import util

from utils_configs.plot.get_columns_from_files import get_columns_from_files
from utils_configs.plot.HEPPlotter import HEPPlotter

import met_ptreg_performance.helpers as helpers

import mc_truth_ptreg_jerc.response_plot.column_utils as cu
from mc_truth_ptreg_jerc.response_plot.write_l2rel_mc_truth import (
    flavour_suffix,
    write_l2relative_txt,
)

DEFAULT_CONFIG = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "mc_truth_configs",
    "mc_truth_config_default.yaml",
)

PUPPI_JET_STRING = r"anti-$k_{T}$ R=0.4 (PUPPI)"


# ---------------------------------------------------------------------------
# command line
# ---------------------------------------------------------------------------
def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Derive the MC truth corrections from the dumped columns",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input-dir",
        default=None,
        help="Directory with the .coffea output of the pocket-coffea run",
    )
    parser.add_argument(
        "-o", "--output", default="./mc_truth_plots", help="Output directory"
    )
    parser.add_argument(
        "-c",
        "--config",
        default=None,
        help="YAML configuration; merged on top of mc_truth_config_default.yaml",
    )
    parser.add_argument(
        "-w", "--workers", type=int, default=1, help="Workers used to draw the plots"
    )
    parser.add_argument(
        "-m",
        "--max-parquet",
        type=int,
        default=None,
        help="Maximum number of parquet files read per dataset (quick tests)",
    )
    parser.add_argument(
        "-l",
        "--load",
        default=None,
        help="Reload a previously produced results file instead of the columns",
    )
    parser.add_argument(
        "--histo",
        action="store_true",
        help="Also plot the 1D response histogram of every bin",
    )
    parser.add_argument(
        "--max-histo-plots",
        type=int,
        default=2000,
        help="Safety limit on the number of response histogram plots",
    )
    parser.add_argument(
        "--no-save-histograms",
        action="store_true",
        help=(
            "Do not store the response histograms in the results file "
            "(smaller output, but --load cannot redraw them)"
        ),
    )
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Only derive the corrections (fits, JSON and txt files), no plots",
    )
    parser.add_argument(
        "--no-flavour-split",
        action="store_true",
        help="Only process the inclusive flavour",
    )
    parser.add_argument(
        "--novars",
        action="store_true",
        help="Old save format, without the variations stored in the coffea files",
    )
    parser.add_argument(
        "--test", action="store_true", help="Use the reduced test binning"
    )
    return parser.parse_args(argv)


# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------
def get_style(cfg, response_var, key, default=None):
    """Look up color/marker of a response variable by matching its suffix."""
    best = None
    for suffix, style in cfg["plot_settings"].items():
        if response_var.endswith(suffix):
            if best is None or len(suffix) > len(best[0]):
                best = (suffix, style)
    return best[1].get(key, default) if best else default


def flavour_style(cfg, flavour, key, default=None):
    return cfg["flavour_split"].get(f"{key}s", {}).get(flavour, default)


def flavour_label(cfg, flavour):
    return cfg["flavour_split"].get("labels", {}).get(flavour, flavour)


def axis_of_kind(cfg, bin_vars, kind):
    """Return the first bin variable of a given kind (e.g. the flavour axis)."""
    for bin_var in bin_vars:
        if cfg["bin_variables"][bin_var].get("kind") == kind:
            return bin_var
    return None


def x_axis_of(cfg, bin_vars):
    for bin_var in bin_vars:
        if cfg["bin_variables"][bin_var].get("x_variable", False):
            return bin_var
    raise ValueError(f"No x_variable among the bin variables {bin_vars}")


def eta_axis_of(cfg, bin_vars):
    for bin_var in bin_vars:
        binning = cfg["bin_variables"][bin_var]
        if binning.get("kind") != "flavour" and not binning.get("x_variable", False):
            return bin_var
    raise ValueError(f"No eta-like axis among the bin variables {bin_vars}")


def select_flavour(h, flavour_axis, flavour, flavour_groups):
    """Project a histogram onto one flavour (or sum them all for 'inclusive')."""
    if flavour_axis is None:
        return h
    if flavour == "inclusive":
        return h[{flavour_axis: sum}]
    return h[{flavour_axis: flavour_groups.index(flavour)}]


def run_plotter(plotter):
    plotter.run()


def run_plotters(plotters, workers, what):
    if not plotters:
        return
    log.info("Drawing %d %s plots...", len(plotters), what)
    if workers == 1:
        for plotter in plotters:
            plotter.run()
    else:
        with Pool(workers) as pool:
            pool.map(run_plotter, plotters)


# ---------------------------------------------------------------------------
# per-bin computation
# ---------------------------------------------------------------------------
def compute_binned_results(cfg, response_var, h_response, h_mean_x, flavours,
                           flavour_groups):
    """Median, inverse median, resolution and mean reco pT for every bin.

    Returns ``{flavour: {quantity: array of shape (n_eta, n_pt)}}``.
    """
    var_cfg = cfg["response_variables"][response_var]
    bin_vars = var_cfg["bin_vars"]
    flavour_axis = axis_of_kind(cfg, bin_vars, "flavour")
    eta_axis = eta_axis_of(cfg, bin_vars)
    pt_axis = x_axis_of(cfg, bin_vars)

    eta_edges = np.asarray(cfg["bin_variables"][eta_axis]["bin_edges"], dtype=float)
    n_eta = len(eta_edges) - 1
    n_pt = len(cfg["bin_variables"][pt_axis]["bin_edges"]) - 1

    eta_max = var_cfg.get("eta_max")
    min_events = cfg.get("min_events_per_bin", 20)
    median_min = cfg.get("median_min", 0.0)
    estimator = cfg.get("resolution_estimator", "quantile")
    conf_level = cfg.get("ci_conf_level", 0.87)

    results = {}
    for flavour in flavours:
        h_flavour = select_flavour(h_response, flavour_axis, flavour, flavour_groups)
        # bring the axes in the (eta, pt, response) order expected below
        h_flavour = h_flavour.project(eta_axis, pt_axis, response_var)
        counts_view = h_flavour.view()
        response_edges = h_flavour.axes[-1].edges

        h_mean_flavour = None
        if h_mean_x is not None:
            h_mean_flavour = select_flavour(
                h_mean_x, flavour_axis, flavour, flavour_groups
            ).project(eta_axis, pt_axis)

        arrays = {
            key: np.full((n_eta, n_pt), np.nan)
            for key in (
                "median",
                "median_err",
                "inv_median",
                "inv_median_err",
                "resolution",
                "resolution_err",
                "x_mean",
                "x_mean_err",
            )
        }
        arrays["n_events"] = np.zeros((n_eta, n_pt))

        for i_eta in range(n_eta):
            if eta_max is not None and min(
                abs(eta_edges[i_eta]), abs(eta_edges[i_eta + 1])
            ) >= eta_max:
                continue
            for i_pt in range(n_pt):
                counts = np.asarray(counts_view[i_eta, i_pt, :], dtype=float)
                n_entries = counts.sum()
                arrays["n_events"][i_eta, i_pt] = n_entries
                if n_entries < min_events:
                    continue

                median, median_err = cu.median_and_error(counts, response_edges)
                if not np.isfinite(median) or median <= median_min:
                    continue

                arrays["median"][i_eta, i_pt] = median
                arrays["median_err"][i_eta, i_pt] = median_err
                arrays["inv_median"][i_eta, i_pt] = 1.0 / median
                arrays["inv_median_err"][i_eta, i_pt] = median_err / median**2

                if estimator == "confidence":
                    width, width_err = cu.confidence_width(
                        counts,
                        (response_edges[:-1] + response_edges[1:]) / 2,
                        conf_level,
                    )
                else:
                    width, width_err = cu.quantile_width(counts, response_edges)
                if cfg.get("normalize_resolution_by_median", True):
                    width, width_err = width / median, width_err / median
                arrays["resolution"][i_eta, i_pt] = width
                arrays["resolution_err"][i_eta, i_pt] = width_err

                if h_mean_flavour is not None:
                    mean, mean_err = cu.mean_from_mean_storage(
                        h_mean_flavour, (i_eta, i_pt)
                    )
                    arrays["x_mean"][i_eta, i_pt] = mean
                    arrays["x_mean_err"][i_eta, i_pt] = mean_err

        results[flavour] = arrays
        log.info(
            "  %s [%s]: %d bins with a valid median",
            response_var,
            flavour,
            int(np.sum(np.isfinite(arrays["median"]))),
        )
    return results


def fit_inverse_medians(cfg, results, flavours):
    """Fit the inverse median versus the mean reco pT for every (variable, flavour, eta)."""
    max_order = cfg["num_params"] - 3
    fit_results = {}

    for response_var, per_flavour in results["per_variable"].items():
        var_cfg = cfg["response_variables"][response_var]
        if not var_cfg.get("derive_correction", False):
            continue
        eta_axis = eta_axis_of(cfg, var_cfg["bin_vars"])
        eta_edges = np.asarray(cfg["bin_variables"][eta_axis]["bin_edges"], dtype=float)

        fit_results[response_var] = {}
        for flavour in flavours:
            if flavour not in per_flavour:
                continue
            arrays = per_flavour[flavour]
            per_eta = {}
            for i_eta in range(len(eta_edges) - 1):
                x = arrays["x_mean"][i_eta]
                y = arrays["inv_median"][i_eta]
                y_err = arrays["inv_median_err"][i_eta]
                mask = (
                    np.isfinite(x)
                    & np.isfinite(y)
                    & np.isfinite(y_err)
                    & (x > cfg.get("fit_pt_min", 1.0))
                )
                if mask.sum() < 2:
                    continue
                fit = cu.fit_pol_log10(
                    x[mask],
                    y[mask],
                    y_err[mask],
                    max_order=max_order,
                    p_value_threshold=cfg.get("fit_p_value_threshold", 0.05),
                )
                if fit is None:
                    log.warning(
                        "Fit failed for %s [%s] eta [%g, %g)",
                        response_var,
                        flavour,
                        eta_edges[i_eta],
                        eta_edges[i_eta + 1],
                    )
                    continue
                per_eta[i_eta] = fit
            fit_results[response_var][flavour] = per_eta
            log.info(
                "  fitted %s [%s]: %d / %d eta bins",
                response_var,
                flavour,
                len(per_eta),
                len(eta_edges) - 1,
            )
    return fit_results


# ---------------------------------------------------------------------------
# outputs
# ---------------------------------------------------------------------------
def to_builtin(obj):
    if isinstance(obj, dict):
        return {str(k): to_builtin(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_builtin(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    return obj


def save_fit_results(fit_results, output_dir, category):
    path = os.path.join(
        output_dir, f"fit_results_inverse_median_{category}.json"
    )
    with open(path, "w") as f:
        json.dump(to_builtin(fit_results), f, indent=2)
    log.info("Fit results saved to %s", path)
    return path


def write_txt_files(cfg, fit_results, year, output_dir, category="baseline"):
    """Write one L2Relative txt file per (correction, flavour)."""
    campaign = cfg["campaign_map"].get(year, year)
    txt_dir = os.path.join(
        output_dir,
        "l2relative_txt" if category == "baseline" else f"l2relative_txt_{category}",
    )
    written = []

    for response_var, per_flavour in fit_results.items():
        var_cfg = cfg["response_variables"][response_var]
        jet_name = var_cfg.get("txt_jet_name")
        if jet_name is None:
            log.warning(
                "No txt_jet_name for '%s': the txt file is not written.", response_var
            )
            continue
        eta_axis = eta_axis_of(cfg, var_cfg["bin_vars"])
        eta_edges = np.asarray(cfg["bin_variables"][eta_axis]["bin_edges"], dtype=float)

        for flavour, per_eta in per_flavour.items():
            file_name = (
                f"{campaign}_{cfg['version']}_MC_L2Relative_"
                f"{jet_name}{flavour_suffix(flavour)}.txt"
            )
            written.append(
                write_l2relative_txt(
                    txt_dir, file_name, eta_edges, per_eta, cfg["num_params"]
                )
            )
    return written


# ---------------------------------------------------------------------------
# plots
# ---------------------------------------------------------------------------
def summary_plotter(graph_data, *, xlabel, ylabel, ratio_label, output_name,
                    annotation, year, y_lim, x_log=True, unity_line=False,
                    x_range=None, annotations=()):
    """Build a HEPPlotter for one summary graph (median / inverse median / ...)."""
    has_ratio = any(
        gd.get("style", {}).get("is_reference", False) for gd in graph_data.values()
    )
    plotter = (
        HEPPlotter()
        .set_plot_config(
            cmstext="Simulation\nPreliminary",
            cmstext_loc=2,
            lumitext=f"{year} (13.6 TeV)" if year else "",
            cmstext_font_size=30,
        )
        .set_labels(xlabel=xlabel, ylabel=ylabel, ratio_label=ratio_label)
        .set_data(graph_data, plot_type="graph")
        .set_options(
            grid=True,
            legend=True,
            legend_loc="upper right",
            x_log=x_log,
            split_legend=False,
            reference_to_den=False,
            ylim_bottom_value=y_lim[0],
            ylim_top_value=y_lim[1],
            **({"set_ylim_ratio": 0.5} if has_ratio else {}),
        )
        .set_output(output_name)
        .add_annotation(
            0.05,
            0.95,
            annotation,
            horizontalalignment="left",
            verticalalignment="top",
        )
    )
    if unity_line and x_range is not None:
        x_line = np.linspace(x_range[0], x_range[1], 10)
        plotter.add_curve(
            x_line, np.ones_like(x_line), color="black", linestyle="--", linewidth=1
        )
    for kwargs in annotations:
        plotter.add_annotation(**kwargs)
    return plotter


def build_graph_entry(x, x_err, y, y_err, *, color, fmt, legend_name, is_reference=False):
    return {
        "data": {
            "x": [np.asarray(x, dtype=float), np.asarray(x_err, dtype=float)],
            "y": [np.asarray(y, dtype=float), np.asarray(y_err, dtype=float)],
        },
        "style": {
            "fmt": fmt,
            "linestyle": None,
            "color": color,
            "markersize": 8,
            "linewidth": 2,
            "legend_name": legend_name,
            "is_reference": is_reference,
        },
    }


def eta_annotation(cfg, eta_edges, i_eta, extra=""):
    text = PUPPI_JET_STRING + "\n" + cu.format_bin_annotation_line(
        eta_edges[i_eta], eta_edges[i_eta + 1], r"$\eta^{\mathrm{reco}}$"
    )
    return text + (f"\n{extra}" if extra else "")


def plot_corrections_per_flavour(cfg, results, fit_results, year, output_dir):
    """One plot per (flavour, eta bin) overlaying all the corrections."""
    plotters = []
    flavours = results["flavours"]

    for flavour in flavours:
        for quantity, ylabel, y_lim, use_reco_x in (
            ("median", "Median response", cfg["y_lim_median"], False),
            (
                "inv_median",
                "1 / Median response",
                cfg["y_lim_inverse_median"],
                True,
            ),
            (
                "resolution",
                (
                    "Response resolution / median"
                    if cfg.get("normalize_resolution_by_median", True)
                    else "Response resolution"
                ),
                cfg["y_lim_resolution"],
                False,
            ),
        ):
            first_var = next(iter(results["per_variable"]))
            eta_axis = eta_axis_of(cfg, cfg["response_variables"][first_var]["bin_vars"])
            eta_edges = np.asarray(
                cfg["bin_variables"][eta_axis]["bin_edges"], dtype=float
            )
            pt_axis = x_axis_of(cfg, cfg["response_variables"][first_var]["bin_vars"])
            pt_edges = np.asarray(cfg["bin_variables"][pt_axis]["bin_edges"], dtype=float)
            pt_centers = (pt_edges[:-1] + pt_edges[1:]) / 2

            for i_eta in range(len(eta_edges) - 1):
                graph_data = {}
                fit_annotations = []
                for response_var, per_flavour in results["per_variable"].items():
                    if flavour not in per_flavour:
                        continue
                    var_cfg = cfg["response_variables"][response_var]
                    arrays = per_flavour[flavour]
                    y = arrays[quantity][i_eta]
                    y_err = arrays[f"{quantity}_err"][i_eta]
                    if not np.any(np.isfinite(y)):
                        continue

                    if use_reco_x:
                        x = arrays["x_mean"][i_eta]
                        x_err = np.zeros_like(x)
                    else:
                        x = pt_centers
                        x_err = np.zeros_like(x)

                    graph_data[response_var] = build_graph_entry(
                        x,
                        x_err,
                        y,
                        y_err,
                        color=get_style(cfg, response_var, "color"),
                        fmt=get_style(cfg, response_var, "fmt", "o"),
                        legend_name=var_cfg.get("legend_name", response_var),
                        is_reference=(
                            quantity == "resolution"
                            and var_cfg.get("is_reference", False)
                        ),
                    )

                    # overlay the polynomial fit of the inverse median
                    fit = (
                        fit_results.get(response_var, {})
                        .get(flavour, {})
                        .get(i_eta)
                    )
                    if use_reco_x and fit is not None:
                        x_fit = np.logspace(
                            np.log10(fit["jet_pt"][0]), np.log10(fit["jet_pt"][1]), 200
                        )
                        graph_data[f"{response_var}_fit"] = {
                            "data": {
                                "x": [x_fit, np.zeros_like(x_fit)],
                                "y": [
                                    cu.evaluate_pol(fit["parameters"], x_fit),
                                    np.zeros(len(x_fit)),
                                ],
                            },
                            "style": {
                                "fmt": "-",
                                "color": get_style(cfg, response_var, "color"),
                                "linewidth": 2,
                                "appear_in_legend": False,
                            },
                        }
                        fit_annotations.append(
                            {
                                "x": 0.05,
                                "y": 0.55 - 0.06 * len(fit_annotations),
                                "text": (
                                    f"pol{fit['pol']}, "
                                    f"$\\chi^2$/ndof = {fit['chi2']:.1f}/{fit['ndof']}"
                                    f", p = {fit['p_value']:.2g}"
                                ),
                                "color": get_style(cfg, response_var, "color"),
                                "horizontalalignment": "left",
                                "verticalalignment": "top",
                                "fontsize": 15,
                            }
                        )

                if not graph_data:
                    continue

                name = (
                    f"{quantity}_{flavour}_jet_eta_"
                    f"{cu.edge_to_string(eta_edges[i_eta])}to"
                    f"{cu.edge_to_string(eta_edges[i_eta + 1])}"
                )
                plotters.append(
                    summary_plotter(
                        graph_data,
                        xlabel=(
                            r"$p_{\mathrm{T}}^{\mathrm{reco}}$ [GeV]"
                            if use_reco_x
                            else r"$p_{\mathrm{T}}^{\mathrm{ptcl}}$ [GeV]"
                        ),
                        ylabel=ylabel,
                        ratio_label="JEC / correction",
                        output_name=os.path.join(output_dir, quantity, name),
                        annotation=eta_annotation(
                            cfg, eta_edges, i_eta, flavour_label(cfg, flavour)
                        ),
                        year=year,
                        y_lim=y_lim,
                        unity_line=(quantity in ("median", "inv_median")),
                        x_range=(pt_edges[0], pt_edges[-1]),
                        annotations=fit_annotations,
                    )
                )
    return plotters


def plot_flavour_comparison(cfg, results, year, output_dir):
    """One plot per (correction, eta bin) overlaying the flavours."""
    plotters = []
    flavours = [f for f in results["flavours"] if f != "inclusive"]
    if not flavours:
        return plotters
    flavours = ["inclusive"] + flavours

    for response_var, per_flavour in results["per_variable"].items():
        var_cfg = cfg["response_variables"][response_var]
        eta_axis = eta_axis_of(cfg, var_cfg["bin_vars"])
        pt_axis = x_axis_of(cfg, var_cfg["bin_vars"])
        eta_edges = np.asarray(cfg["bin_variables"][eta_axis]["bin_edges"], dtype=float)
        pt_edges = np.asarray(cfg["bin_variables"][pt_axis]["bin_edges"], dtype=float)
        pt_centers = (pt_edges[:-1] + pt_edges[1:]) / 2

        for quantity, ylabel, y_lim in (
            ("median", "Median response", cfg["y_lim_median"]),
            (
                "resolution",
                (
                    "Response resolution / median"
                    if cfg.get("normalize_resolution_by_median", True)
                    else "Response resolution"
                ),
                cfg["y_lim_resolution"],
            ),
        ):
            for i_eta in range(len(eta_edges) - 1):
                graph_data = {}
                for flavour in flavours:
                    if flavour not in per_flavour:
                        continue
                    arrays = per_flavour[flavour]
                    y = arrays[quantity][i_eta]
                    y_err = arrays[f"{quantity}_err"][i_eta]
                    if not np.any(np.isfinite(y)):
                        continue
                    graph_data[flavour] = build_graph_entry(
                        pt_centers,
                        np.zeros_like(pt_centers),
                        y,
                        y_err,
                        color=flavour_style(cfg, flavour, "color", "black"),
                        fmt=flavour_style(cfg, flavour, "marker", "o"),
                        legend_name=flavour_label(cfg, flavour),
                        is_reference=(flavour == "inclusive"),
                    )

                if len(graph_data) < 2:
                    continue

                name = (
                    f"{quantity}_flavours_{var_cfg['name_plot']}_jet_eta_"
                    f"{cu.edge_to_string(eta_edges[i_eta])}to"
                    f"{cu.edge_to_string(eta_edges[i_eta + 1])}"
                )
                plotters.append(
                    summary_plotter(
                        graph_data,
                        xlabel=r"$p_{\mathrm{T}}^{\mathrm{ptcl}}$ [GeV]",
                        ylabel=ylabel,
                        ratio_label="all / flavour",
                        output_name=os.path.join(
                            output_dir, f"flavour_{quantity}", name
                        ),
                        annotation=eta_annotation(
                            cfg,
                            eta_edges,
                            i_eta,
                            var_cfg.get("legend_name", response_var),
                        ),
                        year=year,
                        y_lim=y_lim,
                        unity_line=(quantity == "median"),
                        x_range=(pt_edges[0], pt_edges[-1]),
                    )
                )
    return plotters


def plot_response_histograms(cfg, response_h_dict, results, year, output_dir,
                             max_plots):
    """1D response distributions, overlaying the corrections in every bin."""
    plotters = []
    flavours = results["flavours"]
    flavour_groups = results["flavour_groups"]

    # group the response variables that share the same binning
    groups = {}
    for response_var in response_h_dict:
        var_cfg = cfg["response_variables"][response_var]
        key = tuple(
            cfg["bin_variables"][bin_var]["name_plot"]
            for bin_var in var_cfg["bin_vars"]
        )
        groups.setdefault(key, []).append(response_var)

    for group_vars in groups.values():
        reference = group_vars[0]
        ref_cfg = cfg["response_variables"][reference]
        eta_axis = eta_axis_of(cfg, ref_cfg["bin_vars"])
        pt_axis = x_axis_of(cfg, ref_cfg["bin_vars"])
        eta_edges = np.asarray(cfg["bin_variables"][eta_axis]["bin_edges"], dtype=float)
        pt_edges = np.asarray(cfg["bin_variables"][pt_axis]["bin_edges"], dtype=float)

        for flavour in flavours:
            projections = {}
            for response_var in group_vars:
                var_bin_vars = cfg["response_variables"][response_var]["bin_vars"]
                projections[response_var] = select_flavour(
                    response_h_dict[response_var],
                    axis_of_kind(cfg, var_bin_vars, "flavour"),
                    flavour,
                    flavour_groups,
                ).project(
                    eta_axis_of(cfg, var_bin_vars),
                    x_axis_of(cfg, var_bin_vars),
                    response_var,
                )

            for i_eta in range(len(eta_edges) - 1):
                for i_pt in range(len(pt_edges) - 1):
                    if len(plotters) >= max_plots:
                        log.warning(
                            "Reached --max-histo-plots=%d: no more response "
                            "histograms are drawn.",
                            max_plots,
                        )
                        return plotters

                    series = {}
                    for response_var, h_projected in projections.items():
                        var_cfg = cfg["response_variables"][response_var]
                        h_1d = h_projected[i_eta, i_pt, :]
                        if h_1d.values().sum() < cfg.get("min_events_per_bin", 20):
                            continue
                        if cfg.get("rebin_for_plotting", True):
                            h_1d, _ = cu.rebin_histogram(
                                h_1d, max_rebin_factor=cfg.get("rebin_max_factor")
                            )
                        series[response_var] = {
                            "data": h_1d,
                            "style": {
                                "color": get_style(cfg, response_var, "color"),
                                "legend_name": var_cfg.get(
                                    "legend_name", response_var
                                ),
                                "linewidth": 2,
                            },
                        }
                    if not series:
                        continue

                    annotation = (
                        eta_annotation(
                            cfg, eta_edges, i_eta, flavour_label(cfg, flavour)
                        )
                        + "\n"
                        + cu.format_bin_annotation_line(
                            pt_edges[i_pt],
                            pt_edges[i_pt + 1],
                            r"$p_{\mathrm{T}}^{\mathrm{ptcl}}$",
                        )
                    )
                    name = (
                        f"histo_response_{flavour}_jet_eta_"
                        f"{cu.edge_to_string(eta_edges[i_eta])}to"
                        f"{cu.edge_to_string(eta_edges[i_eta + 1])}"
                        f"_jet_gen_pt_{cu.edge_to_string(pt_edges[i_pt])}to"
                        f"{cu.edge_to_string(pt_edges[i_pt + 1])}"
                    )
                    plotters.append(
                        HEPPlotter()
                        .set_plot_config(
                            cmstext="Simulation\nPreliminary",
                            cmstext_loc=2,
                            lumitext=f"{year} (13.6 TeV)" if year else "",
                            cmstext_font_size=30,
                        )
                        .set_data(series, plot_type="1d")
                        .set_labels(
                            xlabel=ref_cfg.get("label", "response"), ylabel="Jets"
                        )
                        .set_output(
                            os.path.join(output_dir, "response_histograms", name)
                        )
                        .set_options(
                            grid=True,
                            legend=True,
                            legend_loc="upper right",
                            split_legend=False,
                        )
                        .add_annotation(
                            0.05,
                            0.95,
                            annotation,
                            horizontalalignment="left",
                            verticalalignment="top",
                        )
                    )
    return plotters


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def flavour_axes(cfg):
    """Names of the bin variables that hold the flavour index."""
    return [
        name
        for name, bin_var_cfg in cfg["bin_variables"].items()
        if bin_var_cfg.get("kind") == "flavour"
    ]


def strip_flavour_axes(cfg):
    """Remove the flavour axes from every variable, keeping only the inclusive sum."""
    to_remove = set(flavour_axes(cfg))
    for section in ("response_variables", "mapping_variables"):
        for var_cfg in cfg.get(section, {}).values():
            var_cfg["bin_vars"] = [
                bin_var
                for bin_var in var_cfg.get("bin_vars", [])
                if bin_var not in to_remove
            ]
    for name in to_remove:
        cfg["bin_variables"].pop(name, None)


def prepare_config(cfg, flavour_groups, no_flavour_split):
    """Fill in the categories of the flavour axes and the flavours to process."""
    enabled = cfg["flavour_split"].get("enabled", False) and not no_flavour_split
    if not enabled:
        strip_flavour_axes(cfg)
        return ["inclusive"]

    for bin_var_cfg in cfg["bin_variables"].values():
        if bin_var_cfg.get("kind") == "flavour":
            bin_var_cfg["categories"] = flavour_groups

    flavours = list(cfg["flavour_split"].get("plot_flavours", ["inclusive"]))
    if "inclusive" not in flavours:
        flavours = ["inclusive"] + flavours
    return flavours


def load_and_fill(cfg, args, flavours, flavour_groups):
    """Read the columns and fill the response and mapping histograms."""
    input_files = [
        os.path.join(args.input_dir, name)
        for name in os.listdir(args.input_dir)
        if name.endswith(".coffea")
    ]
    if not input_files:
        raise SystemExit(f"No .coffea file found in {args.input_dir}")

    log.info("Loading the columns from %d files...", len(input_files))
    cat_col, datasets = get_columns_from_files(
        input_files,
        "nominal",
        None,
        True,
        args.novars,
        max_num_parquet_files=args.max_parquet,
    )

    year = ""
    for dataset in datasets:
        tag = helpers.extract_year_tag(dataset)
        if tag:
            year = ", ".join([year, tag]) if year else tag
    log.info("Datasets: %s (year %s)", datasets, year)

    collections = sorted({name.split("_")[0] for name in flavour_axes(cfg)})

    filled = {}
    for category, columns in cat_col.items():
        log.info("Category '%s': flattening the columns...", category)
        data = cu.flatten_data(columns)
        if collections:
            added = cu.add_flavour_index(data, cfg["flavour_split"], collections)
            if not added:
                log.warning(
                    "No flavour column in the input: falling back to the "
                    "inclusive flavour only. Re-run the analysis with "
                    "partonFlavour/hadronFlavour in the dumped columns to get "
                    "the flavour splitting."
                )
                strip_flavour_axes(cfg)
                collections = []
                flavours[:] = ["inclusive"]

        response_vars = {
            name: var_cfg
            for name, var_cfg in cfg["response_variables"].items()
            if name in data
        }
        if not response_vars:
            log.warning("No response variable found in category '%s'.", category)
            continue

        mapping_vars = {
            name: var_cfg
            for name, var_cfg in cfg["mapping_variables"].items()
            if name in data
        }

        log.info("Filling the response histograms...")
        response_h_dict, _ = cu.create_ND_histo(response_vars, data, cfg["bin_variables"])
        log.info("Filling the mapping histograms...")
        _, mapping_mean_dict = cu.create_ND_histo(
            mapping_vars, data, cfg["bin_variables"]
        )

        filled[category] = {
            "year": year,
            "response_h_dict": response_h_dict,
            "mapping_mean_dict": mapping_mean_dict,
            "response_vars": list(response_vars),
        }
    return filled


def process_category(cfg, args, category, payload, flavours, flavour_groups):
    """Compute, fit, write and plot everything for one category."""
    year = payload["year"]
    response_h_dict = payload["response_h_dict"]
    mapping_mean_dict = payload["mapping_mean_dict"]

    results = {
        "year": year,
        "category": category,
        "flavours": flavours,
        "flavour_groups": flavour_groups,
        "per_variable": {},
    }

    log.info("Computing the per-bin quantities...")
    for response_var in payload["response_vars"]:
        map_x_variable = cfg["response_variables"][response_var].get("map_x_variable")
        h_mean_x = mapping_mean_dict.get(map_x_variable)
        if map_x_variable and h_mean_x is None:
            log.warning(
                "Mapping variable '%s' is missing: the fit of '%s' will use no "
                "reco pT.",
                map_x_variable,
                response_var,
            )
        results["per_variable"][response_var] = compute_binned_results(
            cfg,
            response_var,
            response_h_dict[response_var],
            h_mean_x,
            flavours,
            flavour_groups,
        )

    log.info("Fitting the inverse medians...")
    fit_flavours = [
        flavour
        for flavour in flavours
        if flavour == "inclusive"
        or flavour
        in cfg["flavour_split"].get("derive_correction_flavours", ["inclusive"])
    ]
    results["fit_results"] = fit_inverse_medians(cfg, results, fit_flavours)

    save_fit_results(results["fit_results"], args.output, category)
    if cfg.get("write_txt", True):
        write_txt_files(cfg, results["fit_results"], year, args.output, category)

    output_file = os.path.join(args.output, f"mc_truth_results_{category}.coffea")
    payload_to_save = {"results": results}
    if not args.no_save_histograms:
        # needed to replot the response distributions with --load --histo
        payload_to_save["response_h_dict"] = response_h_dict
    util.save(payload_to_save, output_file)
    log.info("Results saved to %s", output_file)

    return results, response_h_dict


def make_plots(cfg, args, results, response_h_dict):
    if args.no_plots:
        log.info("--no-plots: skipping the plots.")
        return
    year = results["year"]
    plot_dir = os.path.join(args.output, f"plots_{results['category']}")

    plotters = plot_corrections_per_flavour(
        cfg, results, results["fit_results"], year, plot_dir
    )
    run_plotters(plotters, args.workers, "correction summary")

    plotters = plot_flavour_comparison(cfg, results, year, plot_dir)
    run_plotters(plotters, args.workers, "flavour comparison")

    if args.histo and response_h_dict is not None:
        plotters = plot_response_histograms(
            cfg, response_h_dict, results, year, plot_dir, args.max_histo_plots
        )
        run_plotters(plotters, args.workers, "response histogram")


def main(argv=None):
    args = parse_args(argv)
    os.makedirs(args.output, exist_ok=True)
    cu.setup_logging(args.output, "response_mc_truth.log", [log, cu.log])

    cfg = cu.load_yaml_config(
        config_path=args.config,
        default_path=DEFAULT_CONFIG,
        test_override=args.test,
    )
    shutil.copy2(
        args.config or DEFAULT_CONFIG,
        os.path.join(args.output, os.path.basename(args.config or DEFAULT_CONFIG)),
    )

    flavour_groups = cu.flavour_group_names(cfg["flavour_split"])
    flavours = prepare_config(cfg, flavour_groups, args.no_flavour_split)
    log.info("Flavours: %s (groups: %s)", flavours, flavour_groups)

    if args.load:
        log.info("Loading the precomputed results from %s...", args.load)
        loaded = util.load(args.load)
        results = loaded["results"]
        make_plots(cfg, args, results, loaded.get("response_h_dict"))
        log.info("All done!")
        return 0

    if not args.input_dir:
        raise SystemExit("Either -i/--input-dir or -l/--load must be given.")

    filled = load_and_fill(cfg, args, flavours, flavour_groups)
    for category, payload in filled.items():
        log.info("Processing category '%s'...", category)
        results, response_h_dict = process_category(
            cfg, args, category, payload, flavours, flavour_groups
        )
        make_plots(cfg, args, results, response_h_dict)

    log.info("All done!")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
