#!/usr/bin/env python3
"""
Plot JER resolution vs jet eta from a correctionlib JSON file.

This is the JSON counterpart of :mod:`plot_jer_txt`. It loads the
correction from JSON via :mod:`correctionlib`, walks the binning tree to
recover (bin_edges, x_min, x_max, params) rows and then delegates the
actual plotting to :func:`plot_jer_vs_eta_from_data` defined in
``plot_jer_txt.py``.

The JSON is expected to follow the structure produced by
``save_json_resolution`` in ``plot_jer_mc.py``: a single ``Correction``
whose ``data`` is a nested ``Binning`` tree (e.g. JetEta then Rho) with
``Formula`` leaves whose parameters are ``[x_min, x_max, *fit_params]``.

Usage
-----
    python plot_jer_json.py Run3Summer22_V1_NSC_MC_PtResolution_AK4PFPNet.json \
        --pt-values 80 150 300 600 -o plots/

    python plot_jer_json.py /path/to/json/*.json --pt-values 80 150 300 600 -o plots/
"""

import os
import argparse
import logging

import matplotlib
matplotlib.use("Agg")

from correctionlib import schemav2 as cs
from correctionlib.schemav2 import CorrectionSet

from mc_truth_ptreg_jerc.response_plot.plot_jer_txt import (
    plot_jer_vs_eta_from_data,
    _extract_year,
    _extract_jet_algo,
    _YEAR_MAP,
)


log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# correctionlib JSON parsing
# ---------------------------------------------------------------------------

def _strip_clamping_from_formula(expr):
    """
    Inverse of the JERC-json transform in ``save_json_resolution``: replace
    every occurrence of ``max(min(x,[1]),[0])`` with ``x`` and shift parameter
    indices [i] (i>=2) down by 2 so that the resulting expression matches
    the convention used in the JERC txt files (clamping handled externally
    via x_min / x_max).
    """
    import re
    stripped = expr.replace("max(min(x,[1]),[0])", "x")
    stripped = re.sub(
        r"\[(\d+)\]",
        lambda m: f"[{int(m.group(1)) - 2}]",
        stripped,
    )
    return stripped


def _walk_binning(node, bin_vars, current_edges, rows):
    """Recursively walk a Binning tree, appending Formula leaves to *rows*."""
    if isinstance(node, cs.Binning):
        var = node.input
        edges = list(node.edges)
        for i, child in enumerate(node.content):
            lo, hi = float(edges[i]), float(edges[i + 1])
            current_edges[var] = (lo, hi)
            _walk_binning(child, bin_vars, current_edges, rows)
        current_edges.pop(var, None)
        return

    if isinstance(node, cs.Formula):
        params = list(node.parameters or [])
        if len(params) < 2:
            log.error(
                "Formula has fewer than 2 parameters (need x_min, x_max): %s",
                params,
            )
            return
        x_min = float(params[0])
        x_max = float(params[1])
        fit_params = [float(p) for p in params[2:]]
        rows.append(
            dict(
                bin_edges=dict(current_edges),
                x_min=x_min,
                x_max=x_max,
                params=fit_params,
                _expression=node.expression,
                _x_var=node.variables[0] if node.variables else None,
            )
        )
        return

    raise TypeError(
        f"Unsupported node type {type(node).__name__} when walking binning tree. "
        "Only Binning and Formula nodes are handled."
    )


def parse_jer_json(path):
    """
    Parse a correctionlib JSON file containing a single JER resolution
    correction and return (header, rows) in the same shape as
    :func:`plot_jer_txt._parse_jer_txt`.
    """
    with open(path) as fh:
        cset = CorrectionSet.model_validate_json(fh.read())
    if not cset.corrections:
        raise ValueError(f"No corrections found in {path}")

    correction = cset.corrections[0]
    if len(cset.corrections) > 1:
        log.warning(
            "Found %d corrections in %s, using the first one (%s).",
            len(cset.corrections),
            path,
            correction.name,
        )

    bin_vars = []
    node = correction.data
    while isinstance(node, cs.Binning):
        bin_vars.append(node.input)
        node = node.content[0]

    rows = []
    _walk_binning(correction.data, bin_vars, current_edges={}, rows=rows)

    if not rows:
        raise ValueError(f"No Formula leaves found in {path}")

    raw_expr = rows[0]["_expression"]
    x_var = rows[0]["_x_var"]
    if x_var is None:
        raise ValueError(f"Formula leaf has no input variable in {path}")

    formula = _strip_clamping_from_formula(raw_expr)

    # Strip private helper keys before returning
    for r in rows:
        r.pop("_expression", None)
        r.pop("_x_var", None)

    header = dict(bin_vars=bin_vars, x_var=x_var, formula=formula)
    return header, rows


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
        help="JER correctionlib JSON file(s); shell globs like '*.json' accepted",
    )
    parser.add_argument(
        "--pt-values",
        nargs="+", type=float,
        default=[50, 100, 300, 600, 1000, 3000],
        metavar="PT",
        help="JetPt [GeV] values at which to evaluate the resolution "
             "(default: 50 100 300 600 1000 3000)",
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

    for json_path in args.input:
        stem = os.path.splitext(os.path.basename(json_path))[0]
        output_base = os.path.join(args.output, stem)

        header, rows = parse_jer_json(json_path)
        year = _extract_year(json_path, year_map=_YEAR_MAP)
        jet_algo = _extract_jet_algo(json_path)

        plot_jer_vs_eta_from_data(
            header=header,
            rows=rows,
            pt_values=args.pt_values,
            output_base=output_base,
            data_formats=args.formats,
            year=year,
            jet_algo=jet_algo,
        )
        plot_jer_vs_eta_from_data(
            header=header,
            rows=rows,
            pt_values=args.pt_values,
            output_base=output_base + "_abseta_lt_2p5",
            data_formats=args.formats,
            eta_max=2.5,
            year=year,
            jet_algo=jet_algo,
        )


if __name__ == "__main__":
    main()
