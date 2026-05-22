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
    (1.305, "solid"),  # Barrel | Endcap EC1
    (2.5, "solid"),  # EC1 | EC2
    (2.964, "solid"),  # EC2 | Forward HF
    (3.139, "dashed"),  # EC2 | Forward HF
    (5.2, "solid"),  # outer limit
]
# (abs-eta centre, top-line, bottom-line) for region labels
_DETECTOR_REGION_LABELS = [
    (0.65, "Barrel", "BB"),
    (1.9, " ", "EC1"),
    (2.3, "Endcap", ""),
    (2.75, " ", "EC2"),
    (3.0, None, None),
    (4.1, "Forward", "HF"),
]

_YEAR_MAP = {
    "2022_preEE": "Summer22",
    "2022_postEE": "Summer22EE",
    "2023_preBPix": "Summer23",
    "2023_postBPix": "Summer23BPix",
    "2024": "Winter24",
    # Ultra-Legacy — APV variant must come before plain UL16 (longer token wins)
    "UL16_APV": "Summer20UL16APV",
    "UL16": "Summer20UL16",
    "UL17": "Summer20UL17",
    "UL18": "Summer20UL18",
}


PUPPI_JET_STRING = r"anti-$k_{T}$ R=0.4 (PUPPI)"


# Keys are substrings matched against the txt filename stem.
# L2L3Residual must come before L2Residual so the longer key wins.
_JERC_CONFIG = {
    "L2Relative": {
        "ylabel": "Simulated jet response",
        "ylim": (0.8, 1.4),
        "ylim_eta_restricted": (0.8, 1.4),
        "lumitext": None,
        "cmstext": "Simulation Preliminary",
        "fold_eta": True,
        "draw_detector_regions": True,
        "plot_inverse": True,
        "hline": 1.0,
        "plot_vs_pt": False,
    },
    "L2Residual": {
        "ylabel": r"$R^\mathrm{sim}/R^\mathrm{data}$",
        "ylim": (0.7, 1.6),
        "ylim_eta_restricted": (0.7, 1.6),
        "lumitext": None,
        "cmstext": "Preliminary",
        "fold_eta": True,
        "draw_detector_regions": True,
        "plot_inverse": False,
        "hline": 1.0,
        "plot_vs_pt": False,
    },
    "L2L3Residual": {
        "ylabel": "Jet response MC / data",
        "ylim": (0.8, 1.4),
        "ylim_eta_restricted": (0.8, 1.4),
        "lumitext": None,       # None → use data["lumi_text"]
        "cmstext": "Preliminary",
        "fold_eta": True,
        "draw_detector_regions": True,
        "plot_inverse": True,
        "hline": 1.0,
        "plot_vs_pt": True,    # True → also produce a vs-pt smooth-curve plot
        "onlyL3Residual": True,
    },
    "MC_JER": {
        "ylabel": "Jet Energy Resolution",
        "ylim": (0.0, 0.3),
        "ylim_eta_restricted": (0.0, 0.2),
        "lumitext": None,
        "cmstext": "Simulation Preliminary",
        "fold_eta": True,
        "draw_detector_regions": True,
        "plot_inverse": False,
        "hline": None,
        "plot_vs_pt": False,     # JER as f(pT) is natural; enabled by default
    },
    "JERSF": {
        "ylabel": "Jet Energy Resolution Scale Factor",
        "ylim": (0.8, 1.5),
        "ylim_eta_restricted": (0.8, 1.5),
        "lumitext": None,
        "cmstext": "Preliminary",
        "fold_eta": True,
        "draw_detector_regions": True,
        "plot_inverse": False,
        "hline": 1.0,
        "plot_vs_pt": False,
    },
}

_DEFAULT_JERC_CONFIG = {
    "ylabel": "Jet Energy Resolution",
    "ylim": (0.0, 0.3),
    "ylim_eta_restricted": (0.0, 0.2),
    "lumitext": None,
    "cmstext": "Simulation Preliminary",
    "fold_eta": True,
    "draw_detector_regions": True,
    "plot_inverse": False,
    "hline": None,
    "plot_vs_pt": False,
}


def _get_jerc_config(stem):
    """Return the _JERC_CONFIG entry whose key appears as a substring in stem."""
    for key, cfg in _JERC_CONFIG.items():
        if key in stem:
            return cfg
    return _DEFAULT_JERC_CONFIG


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
        n_vals = int(toks[idx])
        idx += 1
        x_min = float(toks[idx])
        idx += 1
        x_max = float(toks[idx])
        idx += 1
        params = [float(v) for v in toks[idx : idx + n_vals - 2]]
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

    return (
        dict(
            bin_vars=bin_vars, x_var=x_var, formula=formula, content_type=content_type
        ),
        rows,
    )


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
        expr = re.sub(rf"\[[a-zA-Z]*{i}\]", repr(float(p)), expr)
    expr = (
        expr.replace("pow(", "np.power(")
        .replace("sqrt(", "np.sqrt(")
        .replace("abs(", "np.abs(")
        .replace("exp(", "np.exp(")
        .replace("log10(", "np.log10(")
        .replace("log(", "np.log(")
    )
    return float(
        eval(expr, {"np": np, "x": float(x), "max": max, "min": min})
    )  # noqa: S307


def _extract_l3residual_formula(formula):
    """Split a compound L2L3Residual formula and return (l3_formula, param_offset).

    Expects formula = L2Relative_factor * L3Residual_factor.  Returns the
    L3Residual factor with its [i] indices re-numbered from 0, plus the
    original index offset so callers can slice per-row params accordingly.
    """
    # Split at top-level '*' (ignore '*' inside parentheses)
    parts = []
    depth = 0
    buf = []
    for ch in formula:
        if ch == "(":
            depth += 1
            buf.append(ch)
        elif ch == ")":
            depth -= 1
            buf.append(ch)
        elif ch == "*" and depth == 0:
            parts.append("".join(buf))
            buf = []
        else:
            buf.append(ch)
    if buf:
        parts.append("".join(buf))

    l3_part = parts[-1]
    indices = [int(m.group(1)) for m in re.finditer(r"\[(\d+)\]", l3_part)]
    if not indices:
        return l3_part, 0
    offset = min(indices)

    def _reindex(m):
        return f"[{int(m.group(1)) - offset}]"

    return re.sub(r"\[(\d+)\]", _reindex, l3_part), offset


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
        "eta": r"$\eta$",
        "rho": r"$\rho$",
        "jetpt": r"Jet $p_\mathrm{T}$",
        "pt": r"$p_\mathrm{T}$",
    }
    return _map.get(var_name.lower(), var_name)


def _extract_year(path, year_map=None):
    """Return a year/era label from the txt file name.

    If *year_map* is provided (dict mapping display-label → season-token,
    e.g. {"2022_preEE": "Summer22"}), the filename is searched for each
    season-token and the matching display-label is returned.
    Falls back to a plain 4-digit year regex when no map is given or matches.

    For DATA files (filename contains "_DATA_"), the CMS run era is also
    extracted from a "Run<YYYY><Era>" token and appended to the label,
    e.g. "Run2023Cv4" → ", Era Cv4".
    """
    basename = os.path.basename(path)
    label = None
    if year_map:
        # Sort by token length descending so "Summer22EE" is tried before "Summer22"
        for lbl, token in sorted(
            year_map.items(), key=lambda kv: len(str(kv[1])), reverse=True
        ):
            if str(token) in basename:
                label = lbl
                break
    if label is None:
        m = re.search(r"(20\d{2})", basename)
        if m:
            label = m.group(1)
    if label is None:
        m = re.search(
            r"(?:Summer|Winter|Fall|Autumn|Spring)(\d{2})", basename, re.IGNORECASE
        )
        if m:
            label = f"20{m.group(1)}"
    if label is not None and "_DATA_" in basename:
        m = re.search(r"Run\d{4}([A-Za-z][A-Za-z0-9]*)", basename)
        if m:
            label = f"{label}, Era {m.group(1)}"
    return label


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


def _combine_file_label(data):
    """Return a short human-readable label for a file in a combined plot.

    Combines the year and formatted jet-algo; falls back to the stem when
    neither can be extracted (e.g. custom txt files with non-standard names).
    """
    year = data.get("year") or ""
    jet_algo = data.get("jet_algo")
    algo_fmt = _format_jet_algo(jet_algo) if jet_algo else ""
    parts = [p for p in [year, algo_fmt] if p]
    return " ".join(parts) if parts else data.get("stem", "")


def _combined_lumi_text(data_list):
    """Return a compact era label for multi-file overlays.

    Collapses to "Run 3 (13.6 TeV)" or "Run 2 (13 TeV)" when all files share
    the same beam energy; falls back to joining individual lumi strings.
    """
    lumi_texts = list(dict.fromkeys(d["lumi_text"] for d in data_list))
    if len(lumi_texts) == 1:
        return lumi_texts[0]
    energies = {"13.6 TeV" if "13.6" in lt else "13 TeV" for lt in lumi_texts}
    if energies == {"13.6 TeV"}:
        return "Run 3 (13.6 TeV)"
    if energies == {"13 TeV"}:
        return "Run 2 (13 TeV)"
    return " / ".join(lumi_texts)


def _make_combined_name(data_list, max_basename=200):
    """Build a filesystem-safe combined output base name (no extension).

    Uses per-file short labels (year + jet algo). Falls back to a short MD5
    hash when the result would still exceed *max_basename* characters.
    """
    import hashlib

    def _slug(d):
        return (
            _combine_file_label(d)
            .replace(", ", "_")
            .replace(" ", "_")
            .replace(",", "")
        )

    labels = [_slug(d) for d in data_list]
    name = "combined_" + "_vs_".join(labels)
    if len(name) <= max_basename:
        return name
    h = hashlib.md5(
        "__vs__".join(d.get("stem", str(i)) for i, d in enumerate(data_list)).encode()
    ).hexdigest()[:8]
    return f"combined_{labels[0]}_and_{len(data_list) - 1}more_{h}"


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
    bin_vars = header["bin_vars"]
    x_var = header["x_var"]
    formula = header["formula"]
    content_type = header.get("content_type", "Resolution")

    cfg = _get_jerc_config(_short_stem(txt_path))
    param_offset = 0
    if cfg.get("onlyL3Residual"):
        formula, param_offset = _extract_l3residual_formula(formula)

    eta_var = next((v for v in bin_vars if "eta" in v.lower()), bin_vars[0])
    # Only treat a variable as rho if its name actually contains "rho";
    # a "JetPhi" second variable must NOT fall into this slot.
    rho_var = next((v for v in bin_vars if "rho" in v.lower()), None)

    year = _extract_year(txt_path, year_map=_YEAR_MAP)
    jet_algo = _extract_jet_algo(txt_path)
    # UL campaigns ran at 13 TeV; Run 3 at 13.6 TeV
    is_ul = "UL" in (year or "")
    e_cm = 13000.0 if is_ul else 13600.0
    lumi_text = f"{year} ({'13' if is_ul else '13.6'} TeV)" if year else "(13.6 TeV)"

    if rho_var:
        rho_bins = sorted({row["bin_edges"][rho_var] for row in rows})
    else:
        rho_bins = [None]
    rho_alphas = np.linspace(0.25, 1.0, len(rho_bins))

    # Accumulate y values keyed by (i_pt, eta_bin, rho_bin) so that any extra
    # bin dimensions (e.g. JetPhi) are averaged out rather than duplicated.
    accum = defaultdict(list)  # key → [y, ...]
    eta_c_of = {}  # eta_bin → eta_c
    i_rho_of = {}  # rho_bin → i_rho
    pt_seen = {}  # i_pt → pt_val

    for i_pt, pt_val in enumerate(pt_values):
        for i_rho, rho_bin in enumerate(rho_bins):
            for row in rows:
                if rho_var and row["bin_edges"].get(rho_var) != rho_bin:
                    continue

                eta_lo, eta_hi = row["bin_edges"][eta_var]
                eta_bin = row["bin_edges"][eta_var]
                eta_c = 0.5 * (eta_lo + eta_hi)
                eta_c_of[eta_bin] = eta_c
                i_rho_of[rho_bin] = i_rho

                # Drop unphysical (pt, eta) combinations: E_jet = pT·cosh(η) < E_beam
                if pt_val * np.cosh(abs(eta_c)) >= e_cm / 2.0:
                    continue

                pt_eval = max(row["x_min"], min(pt_val, row["x_max"]))
                row_params = row["params"][param_offset:]
                try:
                    y = _eval_formula(formula, pt_eval, row_params)
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
        points.append(
            {
                "eta_c": eta_c_of[eta_bin],
                "eta_bin": eta_bin,
                "rho_bin": rho_bin,
                "pt_val": pt_seen[i_pt],
                "i_pt": i_pt,
                "i_rho": i_rho_of[rho_bin],
                "y": y_mean,
            }
        )

    pt_values_used = sorted(pt_seen.items())  # [(i_pt, pt_val), ...]

    return {
        "points": points,
        "rho_bins": rho_bins,
        "rho_alphas": rho_alphas,
        "pt_values_used": pt_values_used,
        "eta_var": eta_var,
        "x_var": x_var,
        "year": year,
        "jet_algo": jet_algo,
        "lumi_text": lumi_text,
        "label": _file_label(txt_path),
        "stem": _short_stem(txt_path),
        "content_type": content_type,
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
        points.append(
            {
                "eta_c": pn["eta_c"],
                "eta_bin": pn["eta_bin"],
                "rho_bin": pn["rho_bin"],
                "pt_val": pn["pt_val"],
                "i_pt": pn["i_pt"],
                "i_rho": pn["i_rho"],
                "y": ratio,
            }
        )
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
        "points": points,
        "rho_bins": data_num["rho_bins"],  # keep numerator's rho structure
        "rho_alphas": data_num["rho_alphas"],
        "pt_values_used": sorted(pt_seen.items()),
        "eta_var": data_num["eta_var"],
        "x_var": data_num["x_var"],
        "year": year_num,
        "jet_algo": data_num["jet_algo"],
        "lumi_text": lumi_text,
        "content_type": data_num.get("content_type", "Resolution"),
    }


# ---------------------------------------------------------------------------
# Layer 2b – combine multiple JERCData dicts into one overlaid dataset
# ---------------------------------------------------------------------------


def combine_jerc_data(data_list):
    """
    Merge multiple JERCData dicts for overlaid plotting on a single canvas.

    Only the central rho bin from each file is retained so the rho dimension
    does not clutter the overlay.  The file index then plays the role that the
    rho-bin index normally plays inside plot_jerc_series:
      - colour  → JetPt value  (unchanged)
      - marker  → file index
      - alpha   → file index (lighter = first file, darker = last file)

    To signal combine-mode to plot_jerc_series, rho_bins is set to a list of
    per-file label strings instead of the usual (lo, hi) tuples.
    """
    n_files = len(data_list)
    file_alphas = np.linspace(0.5, 1.0, n_files) if n_files > 1 else np.array([1.0])
    file_labels = [_combine_file_label(d) for d in data_list]

    all_points = []
    pt_seen = {}
    for i_file, data in enumerate(data_list):
        rho_bins = data["rho_bins"]
        if rho_bins == [None]:
            central_rho_bin = None
        else:
            central_rho_bin = rho_bins[len(rho_bins) // 2]

        for p in data["points"]:
            if rho_bins != [None] and p["rho_bin"] != central_rho_bin:
                continue
            all_points.append({**p, "i_rho": i_file, "rho_bin": file_labels[i_file]})
            pt_seen[p["i_pt"]] = p["pt_val"]

    lumi_text = _combined_lumi_text(data_list)

    algos = list(dict.fromkeys(d["jet_algo"] for d in data_list if d["jet_algo"]))
    jet_algo = algos[0] if len(algos) == 1 else None

    base = data_list[0]
    return {
        "points": all_points,
        "rho_bins": file_labels,  # strings → signals combine-mode
        "rho_alphas": file_alphas,
        "pt_values_used": sorted(pt_seen.items()),
        "eta_var": base["eta_var"],
        "x_var": base["x_var"],
        "year": base["year"],
        "jet_algo": jet_algo,
        "lumi_text": lumi_text,
        "label": " vs ".join(file_labels),
        "stem": "combined",
        "config_stem": base["stem"],  # inherit config (ylim, ylabel…) from first file
        "content_type": base.get("content_type", "Resolution"),
    }


def extract_jerc_vs_pt_curves(txt_path, eta_range, n_pts=300):
    """
    Evaluate the JERC formula on a fine pt grid for a given |eta| window.

    Returns a JERCData-like dict with x_mode="pt" that can be passed
    directly to plot_jerc_series or combine_jerc_vs_pt_curves.

    Parameters
    ----------
    txt_path  : path to the JERC txt file
    eta_range : (eta_lo, eta_hi) — the |eta| window to select
    n_pts     : number of pt grid points (log-spaced 10–5000 GeV)
    """
    eta_lo, eta_hi = eta_range
    stem = _short_stem(txt_path)
    cfg = _get_jerc_config(stem)

    header, rows = _parse_jerc_txt(txt_path)
    bin_vars = header["bin_vars"]
    formula = header["formula"]
    param_offset = 0
    if cfg.get("onlyL3Residual"):
        formula, param_offset = _extract_l3residual_formula(formula)
    eta_var = next((v for v in bin_vars if "eta" in v.lower()), bin_vars[0])
    rho_var = next((v for v in bin_vars if "rho" in v.lower()), None)

    year = _extract_year(txt_path, year_map=_YEAR_MAP)
    jet_algo = _extract_jet_algo(txt_path)
    is_ul = "UL" in (year or "")
    e_cm = 13000.0 if is_ul else 13600.0
    lumi_text = f"{year} ({'13' if is_ul else '13.6'} TeV)" if year else "(13.6 TeV)"

    eta_rows = [
        r for r in rows
        if eta_lo <= abs(0.5 * sum(r["bin_edges"][eta_var])) < eta_hi
    ]
    if not eta_rows:
        log.warning(
            "No eta bins with centre in [%.2f, %.2f) for %s — skipping vs-pt.",
            eta_lo, eta_hi, txt_path,
        )
        return None

    rho_bins_list = sorted({r["bin_edges"][rho_var] for r in eta_rows}) if rho_var else [None]
    rho_alphas = np.linspace(0.25, 1.0, len(rho_bins_list))
    pt_grid = np.logspace(1, np.log10(5000), n_pts)

    points = []
    for i_rho, rho_bin in enumerate(rho_bins_list):
        rho_rows = [
            r for r in eta_rows
            if rho_var is None or r["bin_edges"].get(rho_var) == rho_bin
        ]
        for pt_val in pt_grid:
            y_vals = []
            for row in rho_rows:
                eta_lo_r, eta_hi_r = row["bin_edges"][eta_var]
                eta_c = 0.5 * (eta_lo_r + eta_hi_r)
                if pt_val * np.cosh(abs(eta_c)) >= e_cm / 2.0:
                    continue
                if not (row["x_min"] <= pt_val <= row["x_max"]):
                    continue
                try:
                    y = _eval_formula(formula, pt_val, row["params"][param_offset:])
                except Exception:
                    continue
                if not np.isfinite(y) or y <= 0:
                    continue
                if cfg.get("plot_inverse"):
                    y = 1.0 / y
                y_vals.append(y)
            if y_vals:
                points.append({
                    "x": float(pt_val),
                    "y": float(np.mean(y_vals)),
                    "i_rho": i_rho,
                    "rho_bin": rho_bin,
                    "i_pt": 0,
                    "pt_val": 0.0,
                    "eta_c": 0.0,
                    "eta_bin": (eta_lo, eta_hi),
                })

    return {
        "points": points,
        "rho_bins": rho_bins_list,
        "rho_alphas": rho_alphas,
        "pt_values_used": [],
        "eta_var": eta_var,
        "x_var": "JetPt",
        "year": year,
        "jet_algo": jet_algo,
        "lumi_text": lumi_text,
        "label": stem,
        "stem": stem,
        "content_type": header.get("content_type", "Resolution"),
        "eta_range": eta_range,
        "x_mode": "pt",
    }


def combine_jerc_vs_pt_curves(data_list):
    """
    Merge multiple vs-pt curve dicts for overlaid plotting on one canvas.

    Each file becomes one curve (i_rho = file index), analogous to
    combine_jerc_data for the vs-eta view.
    """
    n_files = len(data_list)
    file_labels = [_combine_file_label(d) for d in data_list]
    # Full opacity: each curve gets its own distinct colour
    file_alphas = np.ones(n_files)

    all_points = []
    for i_file, data in enumerate(data_list):
        for p in data["points"]:
            all_points.append({**p, "i_rho": i_file, "rho_bin": file_labels[i_file]})

    lumi_text = _combined_lumi_text(data_list)

    algos = list(dict.fromkeys(d["jet_algo"] for d in data_list if d["jet_algo"]))
    jet_algo = algos[0] if len(algos) == 1 else None

    base = data_list[0]
    return {
        "points": all_points,
        "rho_bins": file_labels,
        "rho_alphas": file_alphas,
        "pt_values_used": [],
        "eta_var": base["eta_var"],
        "x_var": "JetPt",
        "year": base["year"],
        "jet_algo": jet_algo,
        "lumi_text": lumi_text,
        "label": " vs ".join(file_labels),
        "stem": "combined",
        "config_stem": base["stem"],  # inherit config (ylim, ylabel…) from first file
        "content_type": base.get("content_type", "Resolution"),
        "eta_range": base.get("eta_range"),
        "x_mode": "pt",
    }


# ---------------------------------------------------------------------------
# Layer 2c – eta-range filter (applied after extraction, before plotting)
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
        boundaries = [
            (-eta, ls) for eta, ls in _DETECTOR_REGION_BOUNDARIES
        ] + boundaries

    for eta, ls in boundaries:
        if xmin < eta < xmax:
            ax.axvline(
                eta,
                color="black",
                linestyle=ls,
                linewidth=0.8,
                zorder=1,
                ymin=0.03,
                ymax=0.78,
            )

    # Draw the first boundary at or just beyond the true data edge to mark where the
    # current detector region ends.  Use ax.dataLim (no matplotlib margin) so that a
    # boundary sitting exactly at the data edge is found first rather than being skipped
    # in favour of the next one further out.  Also skip drawing if the boundary already
    # fell inside the padded xlim and was rendered by the main loop above.
    data_xmax = ax.dataLim.x1
    first_outside = next(
        ((eta, ls) for eta, ls in boundaries if eta >= data_xmax and eta > xmin),
        None,
    )
    if first_outside is not None and not (xmin < first_outside[0] < xmax):
        ax.axvline(
            first_outside[0],
            color="black",
            linestyle=first_outside[1],
            linewidth=0.8,
            zorder=1,
            ymin=0.03,
            ymax=0.78,
        )

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
                eta_c,
                0.85,
                text,
                transform=trans,
                ha="center",
                va="top",
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
    lumitext=None,
    cmstext="Simulation Preliminary",
    x_mode="eta",
):
    """
    Turn a JERCData dict into a publication-quality plot and save it.

    Works for both vs-eta (x_mode="eta", default) and vs-pt (x_mode="pt").
    For x_mode="pt" the dict must contain points with an "x" key (pt value);
    use extract_jerc_vs_pt_curves / combine_jerc_vs_pt_curves to build it.

    Parameters
    ----------
    data                  : JERCData dict
    output_base           : output path without extension
    ylabel                : y-axis label string
    annotation_text       : text drawn in the upper-left corner
    ylim                  : (bottom, top) y-axis limits
    data_formats          : iterable of file extensions to save
    eta_max               : if set, draw vertical dashed lines at ±eta_max (eta mode only)
    hline                 : if set, draw a horizontal dashed line at this y value
    fold_eta              : if True, use |eta_c| on the x-axis (eta mode only)
    draw_detector_regions : if True, overlay CMS detector-region lines (eta mode only)
    x_mode                : "eta" (scatter vs eta) or "pt" (smooth curve vs pt)
    """
    points = data["points"]
    rho_bins = data["rho_bins"]
    rho_alphas = data["rho_alphas"]
    pt_values_used = data["pt_values_used"]
    eta_var = data["eta_var"]
    x_var = data["x_var"]
    lumi_text = lumitext if lumitext is not None else data["lumi_text"]

    if not points:
        log.warning("No data points to plot for %s — skipping.", output_base)
        return

    no_rho = rho_bins == [None]

    # Map i_pt → pt_val for label look-ups
    pt_val_of = dict(pt_values_used)

    # Group points and build (x, y, x_err) triples per series key.
    # vs-pt: one smooth curve per rho bin, x = pt value.
    # vs-eta: one scatter series per (i_pt, i_rho), x = eta_c (optionally folded).
    groups = defaultdict(list)
    if x_mode == "pt":
        for p in points:
            groups[(0, p["i_rho"])].append((p["x"], p["y"], 0.0))
    elif fold_eta:
        fold_y = defaultdict(lambda: defaultdict(list))  # (i_pt,i_rho) → |eta_c| → [y]
        fold_xerr = {}  # |eta_c| → x_err
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
        if x_mode == "pt":
            # vs-pt: colour by rho/file index, draw as a smooth line
            base_color = _BASE_COLORS[i_rho % len(_BASE_COLORS)]
            alpha = float(rho_alphas[i_rho]) if i_rho < len(rho_alphas) and not no_rho else 1.0
            color = _rgba(base_color, alpha)
            ser_key = f"rho{i_rho}"
            ser_style = {
                "color": color,
                "fmt": "",
                "markersize": 0,
                "linestyle": "-",
                "linewidth": 1.5,
                "legend_name": "",
                "appear_in_legend": False,
            }
        else:
            # vs-eta: colour by pT, alpha by rho bin
            pt_val = pt_val_of[i_pt]
            base_color = _BASE_COLORS[i_pt % len(_BASE_COLORS)]
            alpha = (
                float(rho_alphas[i_rho]) if i_rho < len(rho_alphas) and not no_rho else 1.0
            )
            color = _rgba(base_color, alpha)
            # When there is no rho binning use a unique marker per pT; otherwise per rho bin
            marker = (
                _MARKERS[i_pt % len(_MARKERS)]
                if no_rho
                else _MARKERS[i_rho % len(_MARKERS)]
            )
            ser_key = f"pt{pt_val:.0f}_rho{i_rho}"
            ser_style = {
                "color": color,
                "fmt": marker,
                "markersize": 5,
                "linestyle": "none",
                "legend_name": f"{x_var} = {pt_val:.0f} GeV",
                "appear_in_legend": False,
            }

        xy_list.sort(key=lambda t: t[0])
        xs = [t[0] for t in xy_list]
        ys = [t[1] for t in xy_list]
        xerrs = [t[2] for t in xy_list]

        series_dict[ser_key] = {
            "data": {"x": [xs, xerrs], "y": ys},
            "style": ser_style,
        }

    if not series_dict:
        log.warning("Empty series for %s — skipping.", output_base)
        return

    if x_mode == "pt":
        xlabel = r"Jet $p_\mathrm{T}$ [GeV]"
    else:
        xlabel = r"$|\eta^{\mathrm{jet}}|$" if fold_eta else _latex_label(eta_var)

    plotter = (
        HEPPlotter("CMS")
        .set_plot_config(
            figsize=None if no_rho else (20, 13),
            cmstext=cmstext,
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
        0.05,
        0.98,
        annotation_text,
        ha="left",
        va="top",
        fontsize="medium",
    )

    fig = plotter.get_figure()
    ax = fig.axes[0]

    if x_mode == "pt":
        ax.set_xscale("log")

    if hline is not None:
        ax.axhline(
            hline, color="black", linestyle="--", linewidth=1.0, alpha=0.6, zorder=0
        )

    if x_mode == "eta" and draw_detector_regions:
        _add_detector_region_lines(ax, symmetric=fold_eta)

    if x_mode == "pt":
        # vs-pt legend: one line per rho bin or file; no pt legend (pt is the x-axis)
        if not no_rho:
            is_combine = any(isinstance(rb, str) for rb in rho_bins)
            rho_handles = [
                mlines.Line2D(
                    [], [],
                    color=_rgba(_BASE_COLORS[i_rho % len(_BASE_COLORS)], float(rho_alphas[i_rho])),
                    linestyle="-",
                    linewidth=2,
                    label=(
                        rho_bin if is_combine
                        else rf"$[{rho_bin[0]:.1f},\ {rho_bin[1]:.1f}]$"
                    ),
                )
                for i_rho, rho_bin in enumerate(rho_bins)
                if rho_bin is not None and i_rho in present_i_rho
            ]
            fig.subplots_adjust(right=0.75)
            ax.legend(
                handles=rho_handles,
                bbox_to_anchor=(1.0, 1.0),
                loc="upper left",
                title="File" if is_combine else r"$\rho$ [GeV/Area]",
                fontsize="x-small",
                title_fontsize="small",
                framealpha=0.85,
            )
    elif no_rho:
        # vs-eta, no rho: single legend with colour+marker per pT value
        pt_handles = [
            mlines.Line2D(
                [],
                [],
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
            bbox_to_anchor=(0.95, 0.0),
            title=None,
            fontsize="x-small",
            framealpha=0.85,
        )
    else:
        # vs-eta with rho: shrink axes and show two external legends
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

        # Detect combine-mode: rho_bins contains label strings instead of (lo,hi) tuples
        is_combine = any(isinstance(rb, str) for rb in rho_bins)

        # Legend part 2: one line/marker per rho bin (or per file in combine-mode)
        rho_handles = [
            mlines.Line2D(
                [],
                [],
                color=(0.3, 0.3, 0.3, float(rho_alphas[i_rho])),
                marker=_MARKERS[i_rho % len(_MARKERS)],
                linestyle="none",
                markersize=8,
                label=(
                    rho_bin
                    if is_combine
                    else rf"$[{rho_bin[0]:.1f},\ {rho_bin[1]:.1f}]$"
                ),
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
            title="File" if is_combine else r"$\rho$ [GeV/Area]",
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
        f"{PUPPI_JET_STRING}\n{_format_jet_algo(jet_algo)}"
        if jet_algo
        else PUPPI_JET_STRING
    )
    cfg = _get_jerc_config(data.get("config_stem", data["stem"]))
    if cfg["plot_inverse"]:
        inv_points = [{**p, "y": 1.0 / p["y"]} for p in data["points"]]
        data = {**data, "points": inv_points}
    ylim = cfg["ylim_eta_restricted"] if eta_max is not None else cfg["ylim"]
    lumi_text = cfg["lumitext"] if cfg["lumitext"] is not None else data["lumi_text"]
    plot_jerc_series(
        data,
        output_base,
        ylabel=cfg["ylabel"],
        annotation_text=annotation,
        ylim=ylim,
        data_formats=data_formats,
        eta_max=eta_max,
        hline=cfg["hline"],
        fold_eta=cfg["fold_eta"],
        draw_detector_regions=cfg["draw_detector_regions"],
        lumitext=lumi_text,
        cmstext=cfg["cmstext"],
    )


def plot_jerc_ratio_vs_eta(
    data_num,
    data_den,
    output_base,
    data_formats=("png", "pdf"),
    eta_max=None,
):
    """Plot the ratio (num / den) for two pre-extracted JERCData dicts."""
    data = compute_ratio_points(
        _filter_eta(data_num, eta_max), _filter_eta(data_den, eta_max)
    )
    if data is None:
        log.warning(
            "No common bins between %s and %s — skipping ratio plot",
            data_num["label"],
            data_den["label"],
        )
        return
    annotation = (
        # f"Numerator:   {data_num['label']}\n"
        # f"Denominator: {data_den['label']}\n"
        f"{PUPPI_JET_STRING}"
    )
    plot_jerc_series(
        data,
        output_base,
        ylabel=f"JER Ratio ({data_num['label']} / {data_den['label']})",
        annotation_text=annotation,
        ylim=(0.5, 1.5),
        data_formats=data_formats,
        eta_max=eta_max,
        hline=1.0,
        fold_eta=True,
        draw_detector_regions="symmetric",
    )


def plot_jerc_vs_pt_data(
    data,
    output_base,
    data_formats=("png", "pdf"),
):
    """Plot a pre-extracted vs-pt JERCData dict (from extract_jerc_vs_pt_curves)."""
    if not data or not data["points"]:
        log.warning("No data points for vs-pt plot %s — skipping.", output_base)
        return
    jet_algo = data["jet_algo"]
    annotation = (
        f"{PUPPI_JET_STRING}\n{_format_jet_algo(jet_algo)}"
        if jet_algo
        else PUPPI_JET_STRING
    )
    eta_range = data.get("eta_range")
    if eta_range is not None:
        eta_lo, eta_hi = eta_range
        eta_label = rf"$|\eta^{{\mathrm{{jet}}}}| \in [{eta_lo:.2g},\,{eta_hi:.2g})$"
        annotation = annotation + "\n" + eta_label
    cfg = _get_jerc_config(data.get("config_stem", data["stem"]))
    lumi_text = cfg["lumitext"] if cfg.get("lumitext") is not None else data["lumi_text"]
    plot_jerc_series(
        data,
        output_base,
        ylabel=cfg["ylabel"],
        annotation_text=annotation,
        ylim=cfg["ylim"],
        data_formats=data_formats,
        hline=cfg["hline"],
        lumitext=lumi_text,
        cmstext=cfg["cmstext"],
        x_mode="pt",
    )


# ---------------------------------------------------------------------------
# vs-pt smooth-curve plot
# ---------------------------------------------------------------------------


def plot_jerc_vs_pt(
    txt_path,
    output_base,
    eta_range,
    data_formats=("png", "pdf"),
    n_pts=300,
):
    """
    Plot JERC quantity vs JetPt as a smooth curve for a selected |eta| range.

    Thin wrapper around extract_jerc_vs_pt_curves + plot_jerc_vs_pt_data.
    """
    data = extract_jerc_vs_pt_curves(txt_path, eta_range, n_pts=n_pts)
    if data is None:
        return
    plot_jerc_vs_pt_data(data, output_base, data_formats=data_formats)


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
        nargs="+",
        type=float,
        default=[15, 30, 90, 300, 1000, 3000],
        metavar="PT",
        help="JetPt [GeV] values at which to evaluate the resolution (default: 15 30 90 300 1000 3000)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=".",
        help="Output directory (default: current directory)",
    )
    parser.add_argument(
        "--formats",
        nargs="+",
        default=["png", "pdf", "svg"],
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
    parser.add_argument(
        "--combine",
        action="store_true",
        help=(
            "Overlay all input files on a single plot.  "
            "Only the central rho bin is used; files are distinguished by marker and alpha.  "
            "Colours still encode JetPt.  Individual per-file plots are still produced."
        ),
    )
    parser.add_argument(
        "--plot-mode",
        choices=["eta", "pt", "both", "auto"],
        default="auto",
        metavar="MODE",
        help=(
            "Which plot type to produce.  "
            "'eta': quantity vs η for different pT values (existing behaviour).  "
            "'pt': smooth curve of quantity vs pT for a selected |η| range (see --pt-eta-range).  "
            "'both': produce both plot types.  "
            "'auto' (default): use the per-correction-type default set in _JERC_CONFIG "
            "(MC_JER produces vs-pT; others produce vs-η)."
        ),
    )
    parser.add_argument(
        "--pt-eta-range",
        nargs=2,
        type=float,
        metavar=("ETA_LO", "ETA_HI"),
        default=[0.0, 1.305],
        help=(
            "Absolute |η| range [lo, hi) used for the vs-pT smooth-curve plot "
            "(default: 0.0 1.305 — central barrel).  "
            "Printed on the plot.  Only relevant when --plot-mode is 'pt', 'both', or 'auto'."
        ),
    )

    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    eta_max_suffix = f"_abseta_lt_{args.eta_max:.2g}".replace(".", "p")
    pt_eta_range = tuple(args.pt_eta_range)
    eta_range_suffix = (
        f"_vspt_eta{args.pt_eta_range[0]:.2g}to{args.pt_eta_range[1]:.2g}"
        .replace(".", "p")
    )

    # Parse and evaluate each file exactly once
    extracted = {
        txt_path: extract_jerc_points(txt_path, args.pt_values)
        for txt_path in args.input
    }

    for txt_path, data in extracted.items():
        output_base = os.path.join(args.output, data["stem"])
        cfg = _get_jerc_config(data["stem"])

        # Determine which plot types to produce for this file
        if args.plot_mode == "eta":
            do_vs_eta, do_vs_pt = True, False
        elif args.plot_mode == "pt":
            do_vs_eta, do_vs_pt = False, True
        elif args.plot_mode == "both":
            do_vs_eta, do_vs_pt = True, True
        else:  # auto: use per-type default from _JERC_CONFIG
            do_vs_pt = cfg.get("plot_vs_pt", False)
            do_vs_eta = not do_vs_pt

        if do_vs_eta:
            plot_jerc_vs_eta(data, output_base, data_formats=args.formats)
            plot_jerc_vs_eta(
                data,
                output_base + eta_max_suffix,
                data_formats=args.formats,
                eta_max=args.eta_max,
            )

        if do_vs_pt:
            plot_jerc_vs_pt(
                txt_path,
                output_base + eta_range_suffix,
                eta_range=pt_eta_range,
                data_formats=args.formats,
            )

    if args.combine:
        if len(args.input) < 2:
            log.warning(
                "--combine requires at least 2 input files; skipping combined plot."
            )
        else:
            data_list = [extracted[p] for p in args.input]
            combined = combine_jerc_data(data_list)
            combined_name = _make_combined_name(data_list)
            output_base = os.path.join(args.output, combined_name)
            plot_jerc_vs_eta(combined, output_base, data_formats=args.formats)
            plot_jerc_vs_eta(
                combined,
                output_base + eta_max_suffix,
                data_formats=args.formats,
                eta_max=args.eta_max,
            )

            # Combined vs-pt plot (respects --plot-mode)
            cfg0 = _get_jerc_config(data_list[0]["stem"])
            do_vs_pt_combined = (
                args.plot_mode in ("pt", "both")
                or (args.plot_mode == "auto" and cfg0.get("plot_vs_pt", False))
            )
            if do_vs_pt_combined:
                vs_pt_curves = [
                    extract_jerc_vs_pt_curves(p, pt_eta_range) for p in args.input
                ]
                vs_pt_curves = [d for d in vs_pt_curves if d is not None]
                if len(vs_pt_curves) >= 2:
                    combined_pt = combine_jerc_vs_pt_curves(vs_pt_curves)
                    plot_jerc_vs_pt_data(
                        combined_pt,
                        output_base + eta_range_suffix,
                        data_formats=args.formats,
                    )

    if not args.no_ratio:
        for path_a, path_b in combinations(args.input, 2):
            data_a, data_b = extracted[path_a], extracted[path_b]
            ratio_base = os.path.join(
                args.output, f"ratio_{data_a['stem']}__over__{data_b['stem']}"
            )
            plot_jerc_ratio_vs_eta(
                data_a, data_b, ratio_base, data_formats=args.formats
            )
            plot_jerc_ratio_vs_eta(
                data_a,
                data_b,
                ratio_base + eta_max_suffix,
                data_formats=args.formats,
                eta_max=args.eta_max,
            )


if __name__ == "__main__":
    main()
