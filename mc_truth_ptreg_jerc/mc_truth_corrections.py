"""Environment-free loading of JME ``L2Relative`` (MC truth) text files.

This module replaces :func:`mc_truth_ptreg_jerc.custom_functions.get_closure_function_information`
for the new column-based MC truth workflow. The old implementation read the
``UPART`` and ``ABS_ETA_INCLUSIVE`` environment variables to decide which eta
rows of the correction file to keep, because the analysis was run one eta
region at a time. The new workflow runs over the full eta range in a single
job, so *all* rows of the file are always kept and no environment variable is
needed.

The returned dictionary has exactly the same structure expected by
:func:`mc_truth_ptreg_jerc.custom_functions.get_mc_truth_corr`, so the two can
be used together unchanged.
"""

import os


def load_l2relative_txt(corr_file):
    """Parse a JME ``L2Relative`` text file into a dictionary.

    Parameters
    ----------
    corr_file : str
        Path to the ``.txt`` correction file.

    Returns
    -------
    dict
        ``function_string``      formula string taken from the file header
        ``corrections_eta_bins`` ``[[eta_low, ...], [eta_high, ...]]``
        ``corrections_phi_bins`` ``[[phi_low, ...], [phi_high, ...]]`` (empty
                                 list when the file is not binned in phi)
        ``num_params``           number of columns declared per row
        ``jet_pt``               ``[[pt_low, ...], [pt_high, ...]]``
        ``params``               list of parameter lists, one per row
    """
    if not os.path.isfile(corr_file):
        raise FileNotFoundError(f"Correction file not found: {corr_file}")

    with open(corr_file, "r") as f:
        lines = f.readlines()

    # The header looks like
    #   {1 JetEta 1 JetPt (<formula>)  Correction L2Relative }
    # optionally with a JetPhi binning in front of JetPt.
    header = lines[0]
    function_string = header.split("JetPt ")[1].split(" Correction")[0].strip()
    has_phi_bins = "JetPhi" in header

    columns = [line.split() for line in lines[1:] if line.strip()]

    corrections_eta_bins = [
        [float(column[0]) for column in columns],
        [float(column[1]) for column in columns],
    ]

    # offset of the remaining fields caused by the optional phi binning
    offset = 2 if has_phi_bins else 0
    corrections_phi_bins = []
    if has_phi_bins:
        corrections_phi_bins = [
            [float(column[2]) for column in columns],
            [float(column[3]) for column in columns],
        ]

    num_params = [int(column[offset + 2]) for column in columns]
    jet_pt = [
        [float(column[offset + 3]) for column in columns],
        [float(column[offset + 4]) for column in columns],
    ]
    params = [
        [float(column[offset + 5 + i]) for i in range(num_params[j] - 2)]
        for j, column in enumerate(columns)
    ]

    return {
        "function_string": function_string,
        "corrections_eta_bins": corrections_eta_bins,
        "corrections_phi_bins": corrections_phi_bins,
        "num_params": num_params,
        "jet_pt": jet_pt,
        "params": params,
    }
