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

_YEAR_MAP = {
    "2022_preEE":    "Summer22",
    "2022_postEE":   "Summer22EE",
    "2023_preBPix":  "Summer23",
    "2023_postBPix": "Summer23BPix",
    "2024":          "Winter24",
}

PUPPI_JET_STRING = r"anti-$k_{T}$ R=0.4 (PUPPI)"


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def _parse_jer_txt(path):
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
        r"\{\s*(\d+)\s+(.*?)\s+\d+\s+(\S+)\s+(.*?)\s+Resolution\}",
        raw[0],
    )
    if not m:
        raise ValueError(f"Cannot parse header: {raw[0]!r}")

    n_bv = int(m.group(1))
    bin_vars = m.group(2).split()[:n_bv]
    x_var = m.group(3)
    formula = m.group(4)

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

    return dict(bin_vars=bin_vars, x_var=x_var, formula=formula), rows


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
    )
    return float(eval(expr, {"np": np, "x": float(x)}))  # noqa: S307


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
    """Compact filesystem-safe label extracted from a txt file path."""
    year = _extract_year(path, year_map=_YEAR_MAP)
    algo = _extract_jet_algo(path)
    parts = [p for p in [year, algo] if p]
    return "_".join(parts) if parts else os.path.splitext(os.path.basename(path))[0]


# ---------------------------------------------------------------------------
# Layer 1 – extract evaluated JER points from a single txt file
# ---------------------------------------------------------------------------

def extract_jer_points(txt_path, pt_values):
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
    JERData dict:
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
    header, rows = _parse_jer_txt(txt_path)
    bin_vars = header["bin_vars"]
    x_var    = header["x_var"]
    formula  = header["formula"]

    eta_var = next((v for v in bin_vars if "eta" in v.lower()), bin_vars[0])
    rho_var = next(
        (v for v in bin_vars if "rho" in v.lower()),
        bin_vars[1] if len(bin_vars) > 1 else None,
    )

    year      = _extract_year(txt_path, year_map=_YEAR_MAP)
    jet_algo  = _extract_jet_algo(txt_path)
    lumi_text = f"{year} (13.6 TeV)" if year else "(13.6 TeV)"

    if rho_var:
        rho_bins = sorted({row["bin_edges"][rho_var] for row in rows})
    else:
        rho_bins = [None]
    rho_alphas = np.linspace(0.25, 1.0, len(rho_bins))

    points = []
    pt_seen = {}  # i_pt → pt_val, only for pt values that produced at least one point

    for i_pt, pt_val in enumerate(pt_values):
        for i_rho, rho_bin in enumerate(rho_bins):
            for row in rows:
                if rho_var and row["bin_edges"].get(rho_var) != rho_bin:
                    continue

                eta_lo, eta_hi = row["bin_edges"][eta_var]
                eta_c = 0.5 * (eta_lo + eta_hi)

                pt_eval = max(row["x_min"], min(pt_val, row["x_max"]))
                try:
                    y = _eval_formula(formula, pt_eval, row["params"])
                except Exception:
                    continue

                if not np.isfinite(y) or y <= 0:
                    continue

                points.append({
                    "eta_c":   eta_c,
                    "eta_bin": row["bin_edges"][eta_var],
                    "rho_bin": rho_bin,
                    "pt_val":  pt_val,
                    "i_pt":    i_pt,
                    "i_rho":   i_rho,
                    "y":       y,
                })
                pt_seen[i_pt] = pt_val

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
    }


# ---------------------------------------------------------------------------
# Layer 2 – compute point-wise ratio between two JERData dicts
# ---------------------------------------------------------------------------

def compute_ratio_points(data_num, data_den):
    """
    Compute the point-wise JER ratio (numerator / denominator) for every
    matching (eta_bin, rho_bin, pt_val) triple that exists in both datasets.

    Parameters
    ----------
    data_num : JERData dict (numerator), from extract_jer_points
    data_den : JERData dict (denominator), from extract_jer_points

    Returns
    -------
    JERData dict in the same format as extract_jer_points, or None if no
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
# Layer 3 – unified plotting function (works for JER or ratio)
# ---------------------------------------------------------------------------

def plot_jer_series(
    data,
    output_base,
    ylabel,
    annotation_text,
    ylim,
    data_formats=("png", "pdf"),
    eta_max=None,
    hline=None,
):
    """
    Turn a JERData dict (from extract_jer_points or compute_ratio_points)
    into a publication-quality plot and save it.

    Parameters
    ----------
    data            : JERData dict
    output_base     : output path without extension
    ylabel          : y-axis label string
    annotation_text : text drawn in the upper-left corner
    ylim            : (bottom, top) y-axis limits
    data_formats    : iterable of file extensions to save
    eta_max         : if set, draw vertical dashed lines at ±eta_max
    hline           : if set, draw a horizontal dashed line at this y value
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

    # Map i_pt → pt_val for label look-ups
    pt_val_of = dict(pt_values_used)

    # Group points by (i_pt, i_rho), sort each group by eta_c
    groups = defaultdict(list)
    for p in points:
        groups[(p["i_pt"], p["i_rho"])].append((p["eta_c"], p["y"]))

    # Only show rho-bin indices that actually appear in the data
    present_i_rho = {p["i_rho"] for p in points}

    # Build HEPPlotter series_dict
    series_dict = {}
    for (i_pt, i_rho), xy_list in groups.items():
        pt_val = pt_val_of[i_pt]
        base_color = _BASE_COLORS[i_pt % len(_BASE_COLORS)]
        alpha      = float(rho_alphas[i_rho]) if i_rho < len(rho_alphas) else 1.0
        color      = _rgba(base_color, alpha)
        marker     = _MARKERS[i_rho % len(_MARKERS)]

        xy_list.sort(key=lambda t: t[0])
        xs = [t[0] for t in xy_list]
        ys = [t[1] for t in xy_list]

        series_dict[f"pt{pt_val:.0f}_rho{i_rho}"] = {
            "data": {"x": xs, "y": ys},
            "style": {
                "color": color,
                "fmt": marker,
                "markersize": 5,
                "legend_name": f"{x_var} = {pt_val:.0f} GeV",
                "appear_in_legend": False,
            },
        }

    if not series_dict:
        log.warning("Empty series for %s — skipping.", output_base)
        return

    plotter = (
        HEPPlotter("CMS")
        .set_plot_config(
            figsize=(20, 13),
            lumitext=lumi_text,
            data_formats=list(data_formats),
        )
        .set_output(output_base)
        .set_labels(
            xlabel=_latex_label(eta_var),
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

    if eta_max is not None:
        for sign in (-1, 1):
            ax.axvline(sign * eta_max, color="gray", linestyle="--", linewidth=1.0, alpha=0.7)

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
# High-level wrappers (accept pre-extracted JERData dicts)
# ---------------------------------------------------------------------------

def plot_jer_vs_eta(
    data,
    output_base,
    data_formats=("png", "pdf"),
    eta_max=None,
):
    """Plot a pre-extracted JERData dict as a resolution-vs-eta plot."""
    data = _filter_eta(data, eta_max)
    if not data["points"]:
        log.warning("No data points after eta filter for %s — skipping.", output_base)
        return
    jet_algo = data["jet_algo"]
    annotation = (
        f"{PUPPI_JET_STRING}\n{_format_jet_algo(jet_algo)}" if jet_algo
        else PUPPI_JET_STRING
    )
    plot_jer_series(
        data, output_base,
        ylabel="Jet Energy Resolution",
        annotation_text=annotation,
        ylim=(0.0, 0.3 if eta_max is None else 0.2),
        data_formats=data_formats,
        eta_max=eta_max,
    )


def plot_jer_ratio_vs_eta(
    data_num,
    data_den,
    output_base,
    data_formats=("png", "pdf"),
    eta_max=None,
):
    """Plot the ratio (num / den) for two pre-extracted JERData dicts."""
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
    plot_jer_series(
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

    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    # Parse and evaluate each file exactly once
    extracted = {
        txt_path: extract_jer_points(txt_path, args.pt_values)
        for txt_path in args.input
    }

    for data in extracted.values():
        output_base = os.path.join(args.output, data["stem"])
        plot_jer_vs_eta(data, output_base, data_formats=args.formats)
        plot_jer_vs_eta(data, output_base + "_abseta_lt_2p5", data_formats=args.formats, eta_max=2.5)

    for path_a, path_b in combinations(args.input, 2):
        data_a, data_b = extracted[path_a], extracted[path_b]
        ratio_base = os.path.join(args.output, f"ratio_{data_a['stem']}__over__{data_b['stem']}")
        plot_jer_ratio_vs_eta(data_a, data_b, ratio_base, data_formats=args.formats)
        plot_jer_ratio_vs_eta(data_a, data_b, ratio_base + "_abseta_lt_2p5", data_formats=args.formats, eta_max=2.5)


if __name__ == "__main__":
    main()
