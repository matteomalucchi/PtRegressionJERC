"""Write JME ``L2Relative`` text files from the MC truth fit results.

The old :mod:`write_l2rel` module re-read the per-eta-bin JSON files that
``response.py`` had just written to disk. Here the fit results are passed
directly as a dictionary, so the txt files are produced in the same pass as the
fits and the flavour splitting comes for free.
"""

import logging
import os

log = logging.getLogger(__name__)


def create_pol_string(num_params):
    """ROOT formula of a polynomial in ``log10(x)`` with ``num_params - 2`` terms."""
    pol_string = "[0]"
    for i in range(1, num_params - 2):
        pol_string += f"+[{i}]*pow(log10(x),{i})"
    return pol_string


def flavour_suffix(flavour):
    """File-name suffix for a flavour ('' for the inclusive sample)."""
    return "" if flavour == "inclusive" else f"_{flavour}Flav"


def write_l2relative_txt(
    output_dir,
    file_name,
    eta_bins,
    fit_results_by_eta,
    num_params,
):
    """Write one L2Relative txt file.

    Parameters
    ----------
    output_dir : str
        Directory where the file is written.
    file_name : str
        Name of the txt file.
    eta_bins : sequence of float
        Eta bin edges; one row is written per bin.
    fit_results_by_eta : dict
        Maps the eta bin index to the fit result dictionary returned by
        :func:`column_utils.fit_pol_log10`. Missing bins are written as a unit
        correction, exactly like the old implementation did.
    num_params : int
        Total number of columns declared per row: the number of polynomial
        parameters plus two (the pT validity range).
    """
    os.makedirs(output_dir, exist_ok=True)
    path = os.path.join(output_dir, file_name)

    n_written, n_missing = 0, 0
    with open(path, "w") as txt_file:
        txt_file.write(
            f"{{1 JetEta 1 JetPt ({create_pol_string(num_params)})"
            "  Correction L2Relative }\n"
        )
        for i in range(len(eta_bins) - 1):
            fit_result = fit_results_by_eta.get(i)
            if not fit_result or not fit_result.get("parameters"):
                # no usable fit: unit correction over a null pT range
                params = ["1"] + ["0"] * (num_params - 3)
                txt_file.write(
                    f" {eta_bins[i]} {eta_bins[i + 1]}    {num_params}  0 0    "
                    + "    ".join(params)
                    + "\n"
                )
                n_missing += 1
                continue

            parameters = list(fit_result["parameters"])
            if len(parameters) > num_params - 2:
                raise ValueError(
                    f"{file_name}: the fit of eta bin {i} has {len(parameters)} "
                    f"parameters but the file declares only {num_params - 2}. "
                    "Increase 'num_params' in the configuration."
                )
            parameters += [0.0] * (num_params - 2 - len(parameters))
            params_string = "".join(f"    {p}" for p in parameters)
            pt_low, pt_high = fit_result["jet_pt"]
            txt_file.write(
                f" {eta_bins[i]} {eta_bins[i + 1]} {num_params}"
                f"    {pt_low}    {pt_high}{params_string}\n"
            )
            n_written += 1

    log.info(
        "Wrote %s (%d fitted eta bins, %d unit-correction bins)",
        path,
        n_written,
        n_missing,
    )
    return path
