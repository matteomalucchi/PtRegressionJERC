import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import os

def plot_nsc_fit_summary(fit_results, bin_var_configs, output_dir, year=""):
    """
    Create summary heatmaps of reduced chi2 and p-value for all NSC fits.

    For each response variable, produces two plots (chi2/ndf and p-value) as
    2D heatmaps when two bin variables are present, or 1D bar charts for one.
    Bins with failed/NaN fits are shown in gray. Red cells flag bad fits
    (chi2/ndf > 3 or p_value < 0.05).
    """

    if not fit_results:
        return

    # Collect all bin variable names that appear in any key (excluding x_var)
    y_bin_var_names = [
        k for k, v in bin_var_configs.items() if not v.get("resolution_x_variable", False)
    ]

    # Group fit_results by response_var
    resp_vars_found = set()
    for key in fit_results:
        for y_var in y_bin_var_names:
            if f"_{y_var}_" in key:
                prefix_end = key.find(f"_{y_var}_")
                resp_vars_found.add(key[:prefix_end])
                break

    for resp_var in sorted(resp_vars_found):
        prefix = f"{resp_var}_"
        matching = {k: v for k, v in fit_results.items() if k.startswith(prefix)}
        if not matching:
            continue

        # Detect which bin variables are present in the keys for this response var
        first_suffix = next(iter(matching.keys()))[len(resp_var):]
        vars_in_key = [k for k in y_bin_var_names if f"_{k}_" in first_suffix]
        if not vars_in_key:
            continue

        # Parse bin edges for each key
        parsed = []
        for key, fit_res in matching.items():
            suffix = key[len(resp_var):]
            edges = {}
            ok = True
            for y_var in vars_in_key:
                marker = f"_{y_var}_"
                idx = suffix.find(marker)
                if idx == -1:
                    ok = False
                    break
                val_start = idx + len(marker)
                val_end = len(suffix)
                for other_var in vars_in_key:
                    other_idx = suffix.find(f"_{other_var}_", val_start)
                    if other_idx != -1 and other_idx < val_end:
                        val_end = other_idx
                try:
                    lo_str, hi_str = suffix[val_start:val_end].split("to", 1)
                    edges[y_var] = (float(lo_str), float(hi_str))
                except ValueError:
                    ok = False
                    break
            if not ok:
                continue
            chi2 = fit_res.get("chi2", np.nan)
            dof = fit_res.get("dof", 0)
            p_val = fit_res.get("p_value", np.nan)
            reduced_chi2 = chi2 / dof if (dof > 0 and np.isfinite(chi2)) else np.nan
            parsed.append((edges, reduced_chi2, p_val))

        if not parsed:
            continue

        def _bin_label(lo, hi):
            def _fmt(v):
                return f"{int(v)}" if float(v) == int(v) else f"{v:.2g}"
            return f"[{_fmt(lo)}, {_fmt(hi)})"

        def _axis_label(var_name):
            return bin_var_configs.get(var_name, {}).get("label", var_name)

        def _fmt_pval(v):
            if not np.isfinite(v):
                return "NaN"
            if v < 0.01:
                s = f"{v:.1e}"
                # strip leading zero from exponent: 5.0e-04 → 5.0e-4
                import re
                s = re.sub(r"e([+-])0*(\d+)$", r"e\1\2", s)
                return s
            return f"{v:.2f}"

        # Continuous p-value colormap on a log scale.
        # LogNorm(1e-6, 1): position = (log10(p) + 6) / 6
        # Anchor positions computed at the key thresholds so the gradient
        # lands exactly where the color transitions should be:
        #   p=1e-6 → 0.000  p=1e-5 → 0.167  p=1e-4 → 0.333
        #   p=1e-3 → 0.500  p=0.01 → 0.667  p=0.05 → 0.783
        #   p=0.50 → 0.950  p=1.00 → 1.000
        _pval_cmap = mcolors.LinearSegmentedColormap.from_list(
            "pval_log",
            [
                (0.000, "#000000"),  # black  (≤ 1e-6, clipped floor)
                (0.167, "#000000"),  # black  (p = 1e-5 — stay black below this)
                (0.333, "#8B0000"),  # dark red    (p = 1e-4)
                (0.500, "#CC0000"),  # red         (p = 1e-3)
                (0.667, "#FF6347"),  # tomato      (p = 0.01)
                (0.783, "#FFD700"),  # gold        (p = 0.05)
                (0.950, "#32CD32"),  # limegreen   (p = 0.5)
                (1.000, "#006400"),  # dark green  (p = 1.0)
            ],
        )
        _pval_cmap.set_bad(color="lightgray")
        _pval_norm = mcolors.LogNorm(vmin=1e-6, vmax=1.0)
        _pval_sm = plt.cm.ScalarMappable(cmap=_pval_cmap, norm=_pval_norm)

        def _text_color_for_bg(rgba):
            """White text on dark backgrounds, black on light — WCAG relative luminance."""
            r, g, b = rgba[:3]
            def _lin(c):
                return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4
            L = 0.2126 * _lin(r) + 0.7152 * _lin(g) + 0.0722 * _lin(b)
            return "white" if L < 0.35 else "black"

        for metric_name, metric_label, vmin, vmax in [
            ("reduced_chi2", r"$\chi^2$/ndf", 0, 5),
            ("p_value", "p-value", 0, 1),
        ]:
            if len(vars_in_key) == 1:
                # 1D horizontal bar chart
                var0 = vars_in_key[0]
                bins0 = sorted(set(e[var0] for e, _, _ in parsed))
                bin_labels = [_bin_label(*b) for b in bins0]

                values = []
                for b in bins0:
                    match = [(r if metric_name == "reduced_chi2" else p)
                             for e, r, p in parsed if e[var0] == b]
                    values.append(match[0] if match else np.nan)

                fig, ax = plt.subplots(figsize=(max(6, len(bins0) * 0.6 + 2), 4))
                _chi2_cmap = plt.get_cmap("RdYlGn_r")
                _chi2_norm = mcolors.Normalize(vmin=0, vmax=5)
                bar_colors = []
                for v in values:
                    if not np.isfinite(v):
                        bar_colors.append("lightgray")
                    elif metric_name == "reduced_chi2":
                        bar_colors.append(_chi2_cmap(_chi2_norm(v)))
                    else:
                        bar_colors.append(_pval_sm.to_rgba(max(v, 1e-7)))

                bars = ax.bar(range(len(bins0)), values, color=bar_colors, edgecolor="black", linewidth=0.5)
                for bar, v in zip(bars, values):
                    label = (_fmt_pval(v) if metric_name == "p_value" else f"{v:.2f}") if np.isfinite(v) else "NaN"
                    bar_top = bar.get_height()
                    ax.text(bar.get_x() + bar.get_width() / 2,
                            bar_top * 1.5 if metric_name == "p_value" else bar_top * 1.01,
                            label, ha="center", va="bottom", fontsize=8)

                if metric_name == "p_value":
                    for thresh, lbl in [(0.05, "p=0.05"), (0.01, "p=0.01"), (1e-4, "p=1e-4"), (1e-5, "p=1e-5")]:
                        ax.axhline(thresh, color=_pval_sm.to_rgba(thresh),
                                   linestyle="--", linewidth=1.2, label=lbl)
                    ax.set_yscale("log")
                    ax.legend(fontsize=9)
                if metric_name == "reduced_chi2":
                    ax.axhline(1, color="black", linestyle="--", linewidth=1)
                    ax.axhline(3, color="red", linestyle="--", linewidth=1)

                ax.set_xticks(range(len(bins0)))
                ax.set_xticklabels(bin_labels, rotation=45, ha="right", fontsize=8)
                ax.set_xlabel(_axis_label(var0), fontsize=10)
                ax.set_ylabel(metric_label, fontsize=10)
                if metric_name == "reduced_chi2":
                    ax.set_ylim(bottom=0, top=vmax)
                else:
                    ax.set_ylim(bottom=1e-6, top=1.5)
                title = f"NSC fit {metric_label} — {resp_var}"
                if year:
                    title += f"  [{year}]"
                ax.set_title(title, fontsize=10)
                fig.tight_layout()

                out_name = f"nsc_fit_summary_{metric_name}_{resp_var}.png"
                fig.savefig(os.path.join(output_dir, out_name), dpi=150)
                plt.close(fig)
                print("Saved NSC fit summary: %s", out_name)

            elif len(vars_in_key) == 2:
                var0, var1 = vars_in_key[0], vars_in_key[1]
                bins0 = sorted(set(e[var0] for e, _, _ in parsed))
                bins1 = sorted(set(e[var1] for e, _, _ in parsed))
                labels0 = [_bin_label(*b) for b in bins0]
                labels1 = [_bin_label(*b) for b in bins1]

                grid = np.full((len(bins1), len(bins0)), np.nan)
                for edges, r_chi2, p_val_v in parsed:
                    i0 = bins0.index(edges[var0])
                    i1 = bins1.index(edges[var1])
                    grid[i1, i0] = r_chi2 if metric_name == "reduced_chi2" else p_val_v

                # Build masked array: NaN → gray overlay
                masked = np.ma.masked_invalid(grid)
                if metric_name == "p_value":
                    cmap = _pval_cmap
                    norm = _pval_norm
                else:
                    cmap = plt.get_cmap("RdYlGn_r").copy()
                    cmap.set_bad(color="lightgray")
                    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

                fig, ax = plt.subplots(figsize=(max(6, len(bins0) * 0.8 + 2),
                                                max(4, len(bins1) * 0.7 + 2)))
                im = ax.imshow(masked, aspect="auto", origin="lower",
                               cmap=cmap, norm=norm,
                               extent=[-0.5, len(bins0) - 0.5, -0.5, len(bins1) - 0.5])

                # Overlay cell values and mark bad fits
                _chi2_sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
                for i1 in range(len(bins1)):
                    for i0 in range(len(bins0)):
                        v = grid[i1, i0]
                        if np.isfinite(v):
                            is_bad = (v > 3) if metric_name == "reduced_chi2" else (v < 0.05)
                            txt = _fmt_pval(v) if metric_name == "p_value" else f"{v:.2f}"
                            cell_color = _pval_sm.to_rgba(max(v, 1e-7)) if metric_name == "p_value" else _chi2_sm.to_rgba(v)
                            ax.text(i0, i1, txt, ha="center", va="center", fontsize=7,
                                    fontweight="bold" if is_bad else "normal",
                                    color=_text_color_for_bg(cell_color))
                            if is_bad:
                                ax.add_patch(plt.Rectangle(
                                    (i0 - 0.5, i1 - 0.5), 1, 1,
                                    fill=False, edgecolor="red", linewidth=2))
                        else:
                            ax.text(i0, i1, "NaN", ha="center", va="center",
                                    fontsize=7, color="dimgray")

                import matplotlib.ticker as mticker
                import re as _re
                cbar = plt.colorbar(im, ax=ax, label=metric_label)
                if metric_name == "p_value":
                    cbar.ax.yaxis.set_major_locator(
                        mticker.LogLocator(base=10, numticks=8)
                    )
                    cbar.ax.yaxis.set_major_formatter(mticker.FuncFormatter(
                        lambda x, _: _re.sub(r"e([+-])0*(\d+)$", r"e\1\2", f"{x:.0e}") if x > 0 else "0"
                    ))
                ax.set_xticks(range(len(bins0)))
                ax.set_xticklabels(labels0, rotation=45, ha="right", fontsize=8)
                ax.set_yticks(range(len(bins1)))
                ax.set_yticklabels(labels1, fontsize=8)
                ax.set_xlabel(_axis_label(var0), fontsize=10)
                ax.set_ylabel(_axis_label(var1), fontsize=10)
                title = f"NSC fit {metric_label} — {resp_var}"
                if year:
                    title += f"  [{year}]"
                ax.set_title(title, fontsize=10)
                fig.tight_layout()

                out_name = f"nsc_fit_summary_{metric_name}_{resp_var}.png"
                fig.savefig(os.path.join(output_dir, out_name), dpi=150)
                plt.close(fig)
                print("Saved NSC fit summary: %s", out_name)

