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
                "Non-finite (NaN/Inf) fit parameter(s) at index %s in line: %r "
                "(bin_edges=%s, x_min=%s, x_max=%s, params=%s) — skipping row",
                nan_param_indices,
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


# ---------------------------------------------------------------------------
# Main plotting function
# ---------------------------------------------------------------------------

def plot_jer_vs_eta(
    txt_path,
    pt_values,
    output_base,
    data_formats=("png", "pdf"),
):
    """
    Read *txt_path* and produce a resolution-vs-eta plot.

    Parameters
    ----------
    txt_path     : path to JERC resolution txt file
    pt_values    : list of JetPt [GeV] at which to evaluate the formula
    output_base  : output path without extension
    data_formats : iterable of file extensions to save
    """
    header, rows = _parse_jer_txt(txt_path)
    bin_vars = header["bin_vars"]
    x_var = header["x_var"]
    formula = header["formula"]

    year = _extract_year(txt_path, year_map=_YEAR_MAP)
    jet_algo = _extract_jet_algo(txt_path)
    lumi_text = f"{year} (13.6 TeV)" if year else "(13.6 TeV)"

    # identify eta and rho variables from bin_vars
    eta_var = next((v for v in bin_vars if "eta" in v.lower()), bin_vars[0])
    rho_var = next(
        (v for v in bin_vars if "rho" in v.lower()),
        bin_vars[1] if len(bin_vars) > 1 else None,
    )

    # sorted unique rho bins (as tuples for hashability)
    if rho_var:
        rho_bins = sorted({row["bin_edges"][rho_var] for row in rows})
    else:
        rho_bins = [None]
    n_rho = len(rho_bins)
    # alpha: low rho → light (0.25), high rho → dark (1.0)
    rho_alphas = np.linspace(0.25, 1.0, n_rho)

    # build HEPPlotter series_dict
    series_dict = {}
    pt_values_used = []  # (original_index, pt_val) only for series that got data

    for i_pt, pt_val in enumerate(pt_values):
        base_color = _BASE_COLORS[i_pt % len(_BASE_COLORS)]
        marker = _MARKERS[i_pt % len(_MARKERS)]
        has_data = False

        for i_rho, rho_bin in enumerate(rho_bins):
            alpha = float(rho_alphas[i_rho])
            color = _rgba(base_color, alpha)

            xs, ys = [], []
            for row in rows:
                # filter by rho bin
                if rho_var and row["bin_edges"].get(rho_var) != rho_bin:
                    continue
                eta_lo, eta_hi = row["bin_edges"][eta_var]
                eta_c = 0.5 * (eta_lo + eta_hi)

                # Clamp pt to the fitted range instead of skipping
                pt_eval = max(row["x_min"], min(pt_val, row["x_max"]))

                try:
                    res = _eval_formula(formula, pt_eval, row["params"])
                except Exception:
                    continue

                if not np.isfinite(res) or res <= 0:
                    continue

                xs.append(eta_c)
                ys.append(res)

            if not xs:
                continue
            has_data = True

            order = np.argsort(xs)
            xs = np.array(xs)[order].tolist()
            ys = np.array(ys)[order].tolist()

            series_key = f"pt{pt_val:.0f}_rho{i_rho}"
            series_dict[series_key] = {
                "data": {"x": xs, "y": ys},
                "style": {
                    "color": color,
                    "fmt": marker,
                    "markersize": 5,
                    "legend_name": f"{x_var} = {pt_val:.0f} GeV",
                    # legend is built manually after get_figure()
                    "appear_in_legend": False,
                },
            }

        if has_data:
            pt_values_used.append((i_pt, pt_val))

    if not series_dict:
        raise RuntimeError(
            "No valid data points were found for the given pt values. "
            "Try different --pt-values within the txt file's JetPt range."
        )

    # ---- HEPPlotter ----
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
            ylabel="Jet Energy Resolution",
        )
        .set_data(series_dict, plot_type="graph")
        .set_options(
            legend=False,           # we add the custom legend below
            set_ylim=True,
            ylim_bottom_value=0.0,
            ylim_top_value=0.3,
            grid=True,
            # enable_watermark=False,
        )
    )

    annotation_text = f"{PUPPI_JET_STRING}\n{_format_jet_algo(jet_algo)}" if jet_algo else PUPPI_JET_STRING
    plotter.add_annotation(
        0.05, 0.98, annotation_text,
        ha="left", va="top",
        fontsize="medium",
    )

    fig = plotter.get_figure()
    ax = fig.axes[0]

    # shrink the axes to leave room for the two external legends on the right
    fig.subplots_adjust(right=0.62)

    # ---- custom two-part legend ----

    # Part 1: one proxy entry per JetPt value (full-alpha base colour)
    pt_handles = []
    for i_pt, pt_val in pt_values_used:
        color_full = _rgba(_BASE_COLORS[i_pt % len(_BASE_COLORS)], 1.0)
        pt_handles.append(
            mlines.Line2D(
                [], [],
                color=color_full,
                marker=_MARKERS[i_pt % len(_MARKERS)],
                linestyle="none",
                markersize=6,
                label=f"{_latex_label(x_var)} = {pt_val:.0f}",
            )
        )

    # Part 2: one proxy entry per Rho bin (greyscale gradient: light → dark)
    rho_handles = []
    for i_rho, rho_bin in enumerate(rho_bins):
        if rho_bin is None:
            continue
        alpha = float(rho_alphas[i_rho])
        rho_handles.append(
            mpatches.Patch(
                facecolor=(0.3, 0.3, 0.3, alpha),
                edgecolor="none",
                label=rf"$\rho \in [{rho_bin[0]:.1f},\ {rho_bin[1]:.1f}]$",
            )
        )

    # Place both legends outside the axes to the right, stacked vertically.
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
        title=r"$\rho$ bins [GeV/Area]" + "\n" + r"(lighter = lower $\rho$)",
        fontsize="x-small",
        title_fontsize="small",
        framealpha=0.85,
    )

    # ---- save ----
    out_dir = os.path.dirname(os.path.abspath(output_base))
    os.makedirs(out_dir, exist_ok=True)
    for fmt in data_formats:
        out_path = f"{output_base}.{fmt}"
        fig.savefig(out_path, bbox_inches="tight", dpi=300, pad_inches=0.05)
        print(f"Saved: {out_path}")
    plt.close(fig)


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
    for txt_path in args.input:
        stem = os.path.splitext(os.path.basename(txt_path))[0]
        output_base = os.path.join(args.output, stem)
        plot_jer_vs_eta(
            txt_path=txt_path,
            pt_values=args.pt_values,
            output_base=output_base,
            data_formats=args.formats,
        )


if __name__ == "__main__":
    main()
