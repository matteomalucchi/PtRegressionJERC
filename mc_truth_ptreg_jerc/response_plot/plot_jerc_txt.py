#!/usr/bin/env python3
"""
Plot JER resolution vs jet eta from a JERC txt file.

The txt file encodes a 2D-binned function (bin vars: JetEta, Rho)
expressed as a function of a third variable (JetPt).

Output:
  - X axis : JetEta bin centre
  - Y axis : resolution evaluated at the chosen JetPt values
  - Colour + marker : JetPt evaluation value
  - Alpha shade : Rho bin  (lighter = lower ρ, darker = higher ρ)

Usage
-----
    # single file
    python plot_jer_txt.py Run3Summer22_V1_NSC_MC_PtResolution_AK4PFPNet.txt \\
        --pt-values 80 150 300 600 -o plots/

    # all txt files in a directory (shell expands the glob)
    python plot_jer_txt.py /path/to/txt/*.txt --pt-values 80 150 300 600 -o plots/
"""

import os
import sys
import re
import argparse
import logging
from collections import defaultdict
from itertools import combinations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import numpy as np

import mplhep as hep
from utils_configs.plot.HEPPlotter import HEPPlotter


log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Style constants
# ---------------------------------------------------------------------------
hep.style.use("CMS")
color_dict = list(hep.style.CMS["axes.prop_cycle"])
_BASE_COLORS = [cycle["color"] for cycle in color_dict]

# _BASE_COLORS = [
#     "#1f77b4",  # blue
#     "#ff7f0e",  # orange
#     "#2ca02c",  # green
#     "#d62728",  # red
#     "#9467bd",  # purple
#     "#8c564b",  # brown
#     "#e377c2",  # pink
# ]
_MARKERS = ["o", "s", "^", "D", "v", "P", "*", "X", "<", ">"]

# CMS detector region boundaries for eta plots
_DETECTOR_REGION_BOUNDARIES = [
    (1.305,  "solid"),   # Barrel | Endcap EC1
    (2.5,  "solid"),   # EC1 | EC2
    (2.964,  "solid"),  # EC2 | Forward HF
    (3.139,  "dashed"),  # EC2 | Forward HF
    (5.2,  "solid"),   # outer limit
]
# (abs-eta centre, top-line, bottom-line) for region labels
_DETECTOR_REGION_LABELS = [
    (0.65, "Barrel",  "BB"),
    (1.9,  " ",  "EC1"),
    (2.3,  "Endcap",  ""),
    (2.75, " ",      "EC2"),
    (3., None,      None),
    (4.1,  "Forward", "HF"),
]

_YEAR_MAP = {
    "2022_preEE":    "Summer22",
    "2022_postEE":   "Summer22EE",
    "2023_preBPix":  "Summer23",
    "2023_postBPix": "Summer23BPix",
    "2024":          "Winter24",
    # Ultra-Legacy — APV variant must come before plain UL16 (longer token wins)
    "UL16_APV":       "Summer20UL16APV",
    "UL16":          "Summer20UL16",
    "UL17":          "Summer20UL17",
    "UL18":          "Summer20UL18",
}


PUPPI_JET_STRING = r"anti-$k_{T}$ R=0.4 (PUPPI)"


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def _parse_jerc_txt(path):
    """
    Parse a JERC-format resolution txt file.

    Header:
        { N_bin_vars BinVar1 [BinVar2 ...] 1 XVar formula Resolution}
    Rows:
        binvar1_lo  binvar1_hi  [...]  n_vals  x_min  x_max  p0  p1  ...

    Returns
    -------
    header : dict  – keys: bin_vars, x_var, formula
    rows   : list of dict – keys: bin_edges, x_min, x_max, params
    """
    with open(path) as fh:
        raw = [ln.strip() for ln in fh if ln.strip()]

    m = re.match(
        r"\{\s*(\d+)\s+(.*?)\s+\d+\s+(\S+)\s+(.*?)\s+(Resolution|Correction\s+L2Relative)\s*\}",
        raw[0],
    )
    if not m:
        raise ValueError(f"Cannot parse header: {raw[0]!r}")

    n_bv = int(m.group(1))
    bin_vars = m.group(2).split()[:n_bv]
    x_var = m.group(3)
    formula = m.group(4)
    content_type = re.sub(r"\s+", " ", m.group(5)).strip()

    rows = []
    for line in raw[1:]:
        toks = line.split()
        idx = 0
        bin_edges = {}
        for var in bin_vars:
            bin_edges[var] = (float(toks[idx]), float(toks[idx + 1]))
            idx += 2
        n_vals = int(toks[idx]); idx += 1
        x_min = float(toks[idx]); idx += 1
        x_max = float(toks[idx]); idx += 1
        params = [float(v) for v in toks[idx: idx + n_vals - 2]]
        nan_param_indices = [i for i, p in enumerate(params) if not np.isfinite(p)]
        if nan_param_indices:
            log.error(
                "Non-finite (NaN/Inf) fit parameter(s) at index %s in file '%s' in line: %r "
                "(bin_edges=%s, x_min=%s, x_max=%s, params=%s) — skipping row",
                nan_param_indices,
                path,
                line,
                bin_edges,
                x_min,
                x_max,
                params,
            )
            continue
        rows.append(dict(bin_edges=bin_edges, x_min=x_min, x_max=x_max, params=params))

    return dict(bin_vars=bin_vars, x_var=x_var, formula=formula, content_type=content_type), rows


# ---------------------------------------------------------------------------
# Formula evaluation
# ---------------------------------------------------------------------------

def _eval_formula(formula_str, x, params):
    """
    Evaluate a ROOT TFormula string at scalar x with the given params list.

    Substitutes [i] → params[i] and maps ROOT math to numpy.
    """
    expr = formula_str
    for i, p in enumerate(params):
        expr = expr.replace(f"[{i}]", repr(float(p)))
    expr = (
        expr.replace("pow(", "np.power(")
           .replace("sqrt(", "np.sqrt(")
           .replace("abs(", "np.abs(")
           .replace("exp(", "np.exp(")
           .replace("log10(", "np.log10(")
           .replace("log(", "np.log(")
    )
    return float(eval(expr, {"np": np, "x": float(x), "max": max, "min": min}))  # noqa: S307


# ---------------------------------------------------------------------------
# Colour helper
# ---------------------------------------------------------------------------

def _rgba(base, alpha):
    """Return (r, g, b, alpha) from a colour spec and an alpha in [0, 1]."""
    r, g, b, _ = mcolors.to_rgba(base)
    return (r, g, b, float(alpha))


# ---------------------------------------------------------------------------
# Label helpers
# ---------------------------------------------------------------------------

def _latex_label(var_name):
    """Map a JERC variable name to a LaTeX-formatted axis/legend label."""
    _map = {
        "jeteta": r"Jet $\eta$",
        "eta":    r"$\eta$",
        "rho":    r"$\rho$",
        "jetpt":  r"Jet $p_\mathrm{T}$",
        "pt":     r"$p_\mathrm{T}$",
    }
    return _map.get(var_name.lower(), var_name)


def _extract_year(path, year_map=None):
    """Return a year/era label from the txt file name.

    If *year_map* is provided (dict mapping display-label → season-token,
    e.g. {"2022_preEE": "Summer22"}), the filename is searched for each
    season-token and the matching display-label is returned.
    Falls back to a plain 4-digit year regex when no map is given or matches.
    """
    basename = os.path.basename(path)
    if year_map:
        # Sort by token length descending so "Summer22EE" is tried before "Summer22"
        for label, token in sorted(year_map.items(), key=lambda kv: len(str(kv[1])), reverse=True):
            if str(token) in basename:
                return label
    m = re.search(r"(20\d{2})", basename)
    if m:
        return m.group(1)
    m = re.search(r"(?:Summer|Winter|Fall|Autumn|Spring)(\d{2})", basename, re.IGNORECASE)
    if m:
        return f"20{m.group(1)}"
    return None


def _extract_jet_algo(path):
    """Return the jet algorithm name (e.g. AK4PFPuppi) from the txt file name, or None."""
    basename = os.path.splitext(os.path.basename(path))[0]
    m = re.search(r"(AK\d+\w+)", basename, re.IGNORECASE)
    return m.group(1) if m else None


def _format_jet_algo(algo):
    """Format a raw jet algorithm string for display.

    Strips the 'AK4PF' prefix, then applies readable substitutions:
      'Plus'     -> ' incl. '
      'Neutrino' -> 'neutrino'
    """
    s = re.sub(r"^AK\d+PF", "", algo)
    s = s.replace("Plus", " incl. ")
    s = s.replace("Neutrino", "neutrino")
    return s


def _file_label(path):
    """Return 'algo-fmt' label for use in ratio plot titles and filenames."""
    algo = _extract_jet_algo(path)
    algo_fmt = _format_jet_algo(algo) if algo else None
    return algo_fmt if algo_fmt else os.path.splitext(os.path.basename(path))[0]


def _short_stem(path):
    """Return the input filename stem (no extension) as the output base name."""
    return os.path.splitext(os.path.basename(path))[0]


# ---------------------------------------------------------------------------
# Layer 1 – extract evaluated JERC points from a single txt file
# ---------------------------------------------------------------------------

def extract_jerc_points(txt_path, pt_values):
    """
    Parse a JERC txt file and evaluate the resolution formula at every
    (JetPt, eta-bin, rho-bin) combination.  All eta bins are included;
    use _filter_eta() afterwards to restrict the range for a specific plot.

    Parameters
    ----------
    txt_path  : path to a JERC resolution txt file
    pt_values : list of nominal JetPt [GeV] values to evaluate

    Returns
    -------
    JERCData dict:
      "points"         – list of point dicts, each with keys:
                           eta_c, eta_bin, rho_bin, pt_val, i_pt, i_rho, y
      "rho_bins"       – sorted list of unique rho-bin tuples (or [None])
      "rho_alphas"     – np.array of alpha values, one per rho bin
      "pt_values_used" – [(i_pt, pt_val), ...] sorted by i_pt
      "eta_var"        – name of the eta bin variable
      "x_var"          – name of the x (JetPt) variable
      "year"           – year/era string extracted from the filename, or None
      "jet_algo"       – raw jet-algorithm string (e.g. "AK4PFPNet"), or None
      "lumi_text"      – ready-to-use lumi string for HEPPlotter
      "label"          – human-readable file label (year + algo) for annotations
      "stem"           – filesystem-safe label for output filenames
    """
    header, rows = _parse_jerc_txt(txt_path)
    bin_vars     = header["bin_vars"]
    x_var        = header["x_var"]
    formula      = header["formula"]
    content_type = header.get("content_type", "Resolution")

    eta_var = next((v for v in bin_vars if "eta" in v.lower()), bin_vars[0])
    # Only treat a variable as rho if its name actually contains "rho";
    # a "JetPhi" second variable must NOT fall into this slot.
    rho_var = next((v for v in bin_vars if "rho" in v.lower()), None)

    year      = _extract_year(txt_path, year_map=_YEAR_MAP)
    jet_algo  = _extract_jet_algo(txt_path)
    # UL campaigns ran at 13 TeV; Run 3 at 13.6 TeV
    is_ul     = "UL" in (year or "")
    e_cm      = 13000.0 if is_ul else 13600.0
    lumi_text = f"{year} ({'13' if is_ul else '13.6'} TeV)" if year else "(13.6 TeV)"

    if rho_var:
        rho_bins = sorted({row["bin_edges"][rho_var] for row in rows})
    else:
        rho_bins = [None]
    rho_alphas = np.linspace(0.25, 1.0, len(rho_bins))

    # Accumulate y values keyed by (i_pt, eta_bin, rho_bin) so that any extra
    # bin dimensions (e.g. JetPhi) are averaged out rather than duplicated.
    accum    = defaultdict(list)   # key → [y, ...]
    eta_c_of = {}                  # eta_bin → eta_c
    i_rho_of = {}                  # rho_bin → i_rho
    pt_seen  = {}                  # i_pt → pt_val

    for i_pt, pt_val in enumerate(pt_values):
        for i_rho, rho_bin in enumerate(rho_bins):
            for row in rows:
                if rho_var and row["bin_edges"].get(rho_var) != rho_bin:
                    continue

                eta_lo, eta_hi = row["bin_edges"][eta_var]
                eta_bin = row["bin_edges"][eta_var]
                eta_c   = 0.5 * (eta_lo + eta_hi)
                eta_c_of[eta_bin] = eta_c
                i_rho_of[rho_bin] = i_rho

                # Drop unphysical (pt, eta) combinations: E_jet = pT·cosh(η) < E_beam
                if pt_val * np.cosh(abs(eta_c)) >= e_cm / 2.0:
                    continue

                pt_eval = max(row["x_min"], min(pt_val, row["x_max"]))
                try:
                    y = _eval_formula(formula, pt_eval, row["params"])
                except Exception:
                    log.error(
                        "Error evaluating formula for pt=%s, eta_bin=%s, rho_bin=%s in file '%s' in line: %r "
                        "(bin_edges=%s, x_min=%s, x_max=%s, params=%s) — skipping point",
                        pt_eval,
                        eta_bin,
                        rho_bin,
                        txt_path,
                        row,
                        row["bin_edges"],
                        row["x_min"],
                        row["x_max"],
                        row["params"],
                    )
                    continue

                if not np.isfinite(y) or y <= 0:
                    continue

                accum[(i_pt, eta_bin, rho_bin)].append(y)
                pt_seen[i_pt] = pt_val

    points = []
    for (i_pt, eta_bin, rho_bin), y_list in accum.items():
        y_mean = float(np.mean(y_list))
        if not np.isfinite(y_mean) or y_mean <= 0:
            continue
        points.append({
            "eta_c":   eta_c_of[eta_bin],
            "eta_bin": eta_bin,
            "rho_bin": rho_bin,
            "pt_val":  pt_seen[i_pt],
            "i_pt":    i_pt,
            "i_rho":   i_rho_of[rho_bin],
            "y":       y_mean,
        })

    pt_values_used = sorted(pt_seen.items())  # [(i_pt, pt_val), ...]

    return {
        "points":          points,
        "rho_bins":        rho_bins,
        "rho_alphas":      rho_alphas,
        "pt_values_used":  pt_values_used,
        "eta_var":         eta_var,
        "x_var":           x_var,
        "year":            year,
        "jet_algo":        jet_algo,
        "lumi_text":       lumi_text,
        "label":           _file_label(txt_path),
        "stem":            _short_stem(txt_path),
        "content_type":    content_type,
    }


# ---------------------------------------------------------------------------
# Layer 2 – compute point-wise ratio between two JERCData dicts
# ---------------------------------------------------------------------------

def compute_ratio_points(data_num, data_den):
    """
    Compute the point-wise JERC ratio (numerator / denominator) for every
    matching (eta_bin, rho_bin, pt_val) triple that exists in both datasets.

    Parameters
    ----------
    data_num : JERCData dict (numerator), from extract_jerc_points
    data_den : JERCData dict (denominator), from extract_jerc_points

    Returns
    -------
    JERCData dict in the same format as extract_jerc_points, or None if no
    common bins are found.  The y values are ratios instead of raw resolutions.
    The colour/marker indices (i_pt, i_rho) are inherited from the numerator
    so that ratio plots use the same visual encoding as the individual plots.
    """
    lut_num = {(p["eta_bin"], p["rho_bin"], p["pt_val"]): p for p in data_num["points"]}
    lut_den = {(p["eta_bin"], p["rho_bin"], p["pt_val"]): p for p in data_den["points"]}

    common_keys = set(lut_num) & set(lut_den)
    if not common_keys:
        log.warning("No common (eta_bin, rho_bin, pt_val) between the two datasets.")
        return None

    points = []
    pt_seen = {}

    for key in common_keys:
        pn = lut_num[key]
        pd = lut_den[key]
        if pd["y"] == 0:
            continue
        ratio = pn["y"] / pd["y"]
        if not np.isfinite(ratio):
            continue
        points.append({
            "eta_c":   pn["eta_c"],
            "eta_bin": pn["eta_bin"],
            "rho_bin": pn["rho_bin"],
            "pt_val":  pn["pt_val"],
            "i_pt":    pn["i_pt"],
            "i_rho":   pn["i_rho"],
            "y":       ratio,
        })
        pt_seen[pn["i_pt"]] = pn["pt_val"]

    if not points:
        return None

    # lumi_text: combine years only when they differ
    year_num = data_num["year"]
    year_den = data_den["year"]
    if year_num == year_den:
        lumi_text = data_num["lumi_text"]
    else:
        lumi_text = f"{year_num or ''}/{year_den or ''} (13.6 TeV)".lstrip("/")

    return {
        "points":          points,
        "rho_bins":        data_num["rho_bins"],   # keep numerator's rho structure
        "rho_alphas":      data_num["rho_alphas"],
        "pt_values_used":  sorted(pt_seen.items()),
        "eta_var":         data_num["eta_var"],
        "x_var":           data_num["x_var"],
        "year":            year_num,
        "jet_algo":        data_num["jet_algo"],
        "lumi_text":       lumi_text,
        "content_type":    data_num.get("content_type", "Resolution"),
    }


# ---------------------------------------------------------------------------
# Layer 2b – eta-range filter (applied after extraction, before plotting)
# ---------------------------------------------------------------------------

def _filter_eta(data, eta_max):
    """Return a shallow copy of *data* with points restricted to |eta_c| < eta_max."""
    if eta_max is None:
        return data
    filtered = [p for p in data["points"] if abs(p["eta_c"]) < eta_max]
    pt_seen = {}
    for p in filtered:
        pt_seen[p["i_pt"]] = p["pt_val"]
    return {**data, "points": filtered, "pt_values_used": sorted(pt_seen.items())}


# ---------------------------------------------------------------------------
# Detector region helpers
# ---------------------------------------------------------------------------

def _add_detector_region_lines(ax, symmetric=False):
    """Draw CMS detector-region vertical lines and text labels on *ax*.

    If *symmetric* is True, mirror lines and labels to negative eta as well.
    """
    import matplotlib.transforms as mtransforms

    xmin, xmax = ax.get_xlim()
    trans = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)

    boundaries = list(_DETECTOR_REGION_BOUNDARIES)
    if symmetric:
        boundaries = [(-eta, ls) for eta, ls in _DETECTOR_REGION_BOUNDARIES] + boundaries

    for eta, ls in boundaries:
        if xmin < eta < xmax:
            ax.axvline(eta, color="black", linestyle=ls, linewidth=0.8, zorder=1,
                       ymin=0.03, ymax=0.78)

    # Draw the first boundary at or just beyond xmax to mark where the current detector region ends
    first_outside = next(
        ((eta, ls) for eta, ls in boundaries if eta >= xmax and eta > xmin),
        None,
    )
    if first_outside is not None:
        ax.axvline(first_outside[0], color="black", linestyle=first_outside[1],
                   linewidth=0.8, zorder=1, ymin=0.03, ymax=0.78)

    label_entries = list(_DETECTOR_REGION_LABELS)
    if symmetric:
        label_entries = [
            (-eta_c, top, bot)
            for eta_c, top, bot in _DETECTOR_REGION_LABELS
            if top is not None or bot is not None
        ] + label_entries

    for eta_c, label_top, label_bot in label_entries:
        if xmin < eta_c < xmax and (label_top is not None or label_bot is not None):
            text = f"{label_top}\n{label_bot}" if label_top else label_bot
            ax.text(
                eta_c, 0.85, text,
                transform=trans,
                ha="center", va="top",
                fontsize=17,
            )


# ---------------------------------------------------------------------------
# Layer 3 – unified plotting function (works for JERC or ratio)
# ---------------------------------------------------------------------------

def plot_jerc_series(
    data,
    output_base,
    ylabel,
    annotation_text,
    ylim,
    data_formats=("png", "pdf"),
    eta_max=None,
    hline=None,
    fold_eta=False,
    draw_detector_regions=False,
):
    """
    Turn a JERCData dict (from extract_jerc_points or compute_ratio_points)
    into a publication-quality plot and save it.

    Parameters
    ----------
    data                  : JERCData dict
    output_base           : output path without extension
    ylabel                : y-axis label string
    annotation_text       : text drawn in the upper-left corner
    ylim                  : (bottom, top) y-axis limits
    data_formats          : iterable of file extensions to save
    eta_max               : if set, draw vertical dashed lines at ±eta_max
    hline                 : if set, draw a horizontal dashed line at this y value
    fold_eta              : if True, use |eta_c| on the x-axis
    draw_detector_regions : if True, overlay CMS detector-region lines and labels
    """
    points          = data["points"]
    rho_bins        = data["rho_bins"]
    rho_alphas      = data["rho_alphas"]
    pt_values_used  = data["pt_values_used"]
    eta_var         = data["eta_var"]
    x_var           = data["x_var"]
    lumi_text       = data["lumi_text"]

    if not points:
        log.warning("No data points to plot for %s — skipping.", output_base)
        return

    no_rho = rho_bins == [None]

    # Map i_pt → pt_val for label look-ups
    pt_val_of = dict(pt_values_used)

    # Group points by (i_pt, i_rho), sort each group by eta_c (or |eta_c|).
    # When folding, average y values of symmetric ±eta bins at the same |eta_c|.
    groups = defaultdict(list)
    if fold_eta:
        fold_y    = defaultdict(lambda: defaultdict(list))  # (i_pt,i_rho) → |eta_c| → [y]
        fold_xerr = {}                                       # |eta_c| → x_err
        for p in points:
            x_val = abs(p["eta_c"])
            eta_lo, eta_hi = p["eta_bin"]
            fold_y[(p["i_pt"], p["i_rho"])][x_val].append(p["y"])
            fold_xerr[x_val] = 0.5 * (eta_hi - eta_lo)
        for key, x_map in fold_y.items():
            for x_val, ys in x_map.items():
                groups[key].append((x_val, float(np.mean(ys)), fold_xerr[x_val]))
    else:
        for p in points:
            eta_lo, eta_hi = p["eta_bin"]
            groups[(p["i_pt"], p["i_rho"])].append(
                (p["eta_c"], p["y"], 0.5 * (eta_hi - eta_lo))
            )

    # Only show rho-bin indices that actually appear in the data
    present_i_rho = {p["i_rho"] for p in points}

    # Build HEPPlotter series_dict
    series_dict = {}
    for (i_pt, i_rho), xy_list in groups.items():
        pt_val = pt_val_of[i_pt]
        base_color = _BASE_COLORS[i_pt % len(_BASE_COLORS)]
        alpha      = float(rho_alphas[i_rho]) if i_rho < len(rho_alphas) and not no_rho else 1.0
        color      = _rgba(base_color, alpha)
        # When there is no rho binning use a unique marker per pT; otherwise per rho bin
        marker     = _MARKERS[i_pt % len(_MARKERS)] if no_rho else _MARKERS[i_rho % len(_MARKERS)]

        xy_list.sort(key=lambda t: t[0])
        xs    = [t[0] for t in xy_list]
        ys    = [t[1] for t in xy_list]
        xerrs = [t[2] for t in xy_list]

        series_dict[f"pt{pt_val:.0f}_rho{i_rho}"] = {
            "data": {"x": [xs, xerrs], "y": ys},
            "style": {
                "color": color,
                "fmt": marker,
                "markersize": 5,
                "linestyle": "none",
                "legend_name": f"{x_var} = {pt_val:.0f} GeV",
                "appear_in_legend": False,
            },
        }

    if not series_dict:
        log.warning("Empty series for %s — skipping.", output_base)
        return

    xlabel = (
        r"$|\eta^{\mathrm{jet}}|$" if fold_eta
        else _latex_label(eta_var)
    )

    plotter = (
        HEPPlotter("CMS")
        .set_plot_config(
            figsize=None if no_rho else (20, 13),
            cmstext="Simulation Preliminary",
            cmstext_font_size=20,
            lumitext=lumi_text,
            lumitext_font_size=20,
            data_formats=list(data_formats),
        )
        .set_output(output_base)
        .set_labels(
            xlabel=xlabel,
            ylabel=ylabel,
        )
        .set_data(series_dict, plot_type="graph")
        .set_options(
            legend=False,
            set_ylim=True,
            ylim_bottom_value=ylim[0],
            ylim_top_value=ylim[1],
            grid=True,
        )
    )

    plotter.add_annotation(
        0.05, 0.98, annotation_text,
        ha="left", va="top",
        fontsize="medium",
    )

    fig = plotter.get_figure()
    ax = fig.axes[0]

    if hline is not None:
        ax.axhline(hline, color="black", linestyle="--", linewidth=1.0, alpha=0.6, zorder=0)

    # if eta_max is not None:
    #     for sign in (-1, 1):
    #         ax.axvline(sign * eta_max, color="gray", linestyle="--", linewidth=1.0, alpha=0.7)

    if draw_detector_regions:
        _add_detector_region_lines(ax, symmetric=(draw_detector_regions == "symmetric"))

    if no_rho:
        # Single legend: colour+marker per pT value, placed inside the plot
        pt_handles = [
            mlines.Line2D(
                [], [],
                color=_rgba(_BASE_COLORS[i_pt % len(_BASE_COLORS)], 1.0),
                marker=_MARKERS[i_pt % len(_MARKERS)],
                linestyle="none",
                markersize=8,
                label=f"$p_{{T}}$ = {pt_val:.0f} GeV",
            )
            for i_pt, pt_val in pt_values_used
        ]
        ax.legend(
            handles=pt_handles,
            loc="lower right",
            title=None,
            fontsize="x-small",
            framealpha=0.85,
        )
    else:
        # shrink axes to leave room for the two external legends on the right
        fig.subplots_adjust(right=0.62)

        # Legend part 1: one colour patch per JetPt value
        pt_handles = [
            mpatches.Patch(
                facecolor=_rgba(_BASE_COLORS[i_pt % len(_BASE_COLORS)], 1.0),
                edgecolor="none",
                label=f"{pt_val:.0f}",
            )
            for i_pt, pt_val in pt_values_used
        ]

        # Legend part 2: one line/marker per rho bin (only those present in the data)
        rho_handles = [
            mlines.Line2D(
                [], [],
                color=(0.3, 0.3, 0.3, float(rho_alphas[i_rho])),
                marker=_MARKERS[i_rho % len(_MARKERS)],
                linestyle="none",
                markersize=8,
                label=rf"$[{rho_bin[0]:.1f},\ {rho_bin[1]:.1f}]$",
            )
            for i_rho, rho_bin in enumerate(rho_bins)
            if rho_bin is not None and i_rho in present_i_rho
        ]

        leg_pt = ax.legend(
            handles=pt_handles,
            bbox_to_anchor=(1.0, 1.0),
            loc="upper left",
            title=f"{_latex_label(x_var)} [GeV]",
            fontsize="x-small",
            title_fontsize="small",
            framealpha=0.85,
        )
        ax.add_artist(leg_pt)

        ax.legend(
            handles=rho_handles,
            bbox_to_anchor=(1.0, 0.0),
            loc="lower left",
            title=r"$\rho$ [GeV/Area]",
            fontsize="x-small",
            title_fontsize="small",
            framealpha=0.85,
        )

    out_dir = os.path.dirname(os.path.abspath(output_base))
    os.makedirs(out_dir, exist_ok=True)
    for fmt in data_formats:
        out_path = f"{output_base}.{fmt}"
        fig.savefig(out_path, bbox_inches="tight", dpi=300, pad_inches=0.05)
        print(f"Saved: {out_path}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# High-level wrappers (accept pre-extracted JERCData dicts)
# ---------------------------------------------------------------------------

def plot_jerc_vs_eta(
    data,
    output_base,
    data_formats=("png", "pdf"),
    eta_max=None,
):
    """Plot a pre-extracted JERCData dict as a resolution-vs-eta (or response-vs-eta) plot."""
    data = _filter_eta(data, eta_max)
    if not data["points"]:
        log.warning("No data points after eta filter for %s — skipping.", output_base)
        return
    jet_algo = data["jet_algo"]
    annotation = (
        f"{PUPPI_JET_STRING}\n{_format_jet_algo(jet_algo)}" if jet_algo
        else PUPPI_JET_STRING
    )
    is_l2 = data.get("content_type", "Resolution") == "Correction L2Relative"
    if is_l2:
        inv_points = [{**p, "y": 1.0 / p["y"]} for p in data["points"]]
        data = {**data, "points": inv_points}
        plot_jerc_series(
            data, output_base,
            ylabel="Simulated jet response",
            annotation_text=annotation,
            ylim=(0.75, 1.25),
            data_formats=data_formats,
            eta_max=eta_max,
            hline=1.0,
            fold_eta=True,
            draw_detector_regions=True,
        )
    else:
        plot_jerc_series(
            data, output_base,
            ylabel="Jet Energy Resolution",
            annotation_text=annotation,
            ylim=(0.0, 0.3 if eta_max is None else 0.2),
            data_formats=data_formats,
            eta_max=eta_max,
            fold_eta=True,
            draw_detector_regions="symmetric",
        )


def plot_jerc_ratio_vs_eta(
    data_num,
    data_den,
    output_base,
    data_formats=("png", "pdf"),
    eta_max=None,
):
    """Plot the ratio (num / den) for two pre-extracted JERCData dicts."""
    data = compute_ratio_points(_filter_eta(data_num, eta_max), _filter_eta(data_den, eta_max))
    if data is None:
        log.warning(
            "No common bins between %s and %s — skipping ratio plot",
            data_num["label"], data_den["label"],
        )
        return
    annotation = (
        # f"Numerator:   {data_num['label']}\n"
        # f"Denominator: {data_den['label']}\n"
        f"{PUPPI_JET_STRING}"
    )
    plot_jerc_series(
        data, output_base,
        ylabel=f"JER Ratio ({data_num['label']} / {data_den['label']})",
        annotation_text=annotation,
        ylim=(0.5, 1.5),
        data_formats=data_formats,
        eta_max=eta_max,
        hline=1.0,
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "input",
        nargs="+",
        help="JERC resolution txt file(s); shell globs like '*.txt' are accepted",
    )
    parser.add_argument(
        "--pt-values",
        nargs="+", type=float,
        default=[50, 100, 300, 600, 1000, 3000,],
        metavar="PT",
        help="JetPt [GeV] values at which to evaluate the resolution (default: 50 100 300 600 1000 3000)",
    )
    parser.add_argument(
        "-o", "--output",
        default=".",
        help="Output directory (default: current directory)",
    )
    parser.add_argument(
        "--formats",
        nargs="+", default=["png", "pdf", "svg"],
        help="Output file formats (default: png pdf svg)",
    )
    parser.add_argument(
        "--no-ratio",
        action="store_true",
        help="Skip ratio plots (useful when only one input file is given)",
    )
    parser.add_argument(
        "--eta-max",
        type=float,
        default=2.5,
        metavar="ETA",
        help="|eta| upper limit for the restricted-range plots (default: 2.5)",
    )

    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    eta_max_suffix = f"_abseta_lt_{args.eta_max:.2g}".replace(".", "p")

    # Parse and evaluate each file exactly once
    extracted = {
        txt_path: extract_jerc_points(txt_path, args.pt_values)
        for txt_path in args.input
    }

    for data in extracted.values():
        output_base = os.path.join(args.output, data["stem"])
        plot_jerc_vs_eta(data, output_base, data_formats=args.formats)
        plot_jerc_vs_eta(data, output_base + eta_max_suffix, data_formats=args.formats, eta_max=args.eta_max)

    if not args.no_ratio:
        for path_a, path_b in combinations(args.input, 2):
            data_a, data_b = extracted[path_a], extracted[path_b]
            ratio_base = os.path.join(args.output, f"ratio_{data_a['stem']}__over__{data_b['stem']}")
            plot_jerc_ratio_vs_eta(data_a, data_b, ratio_base, data_formats=args.formats)
            plot_jerc_ratio_vs_eta(data_a, data_b, ratio_base + eta_max_suffix, data_formats=args.formats, eta_max=args.eta_max)


if __name__ == "__main__":
    main()
