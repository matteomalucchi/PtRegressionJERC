# PtRegressionJERC

Repository with PocketCoffea configurations to compute the MC Truth corrections, JER and MET-Type1 correction for regressed pT jets.

## Setup

### lxplus installation

To setup a local installation on `lxplus`:

```bash
# Clone the fork and checkout the desired branch
git clone https://github.com/PocketCoffea/PocketCoffea.git
cd PocketCoffea

#Enter the Singularity image
apptainer shell --bind /afs -B /cvmfs/cms.cern.ch \
         --bind /tmp  --bind /eos/cms/ -B /etc/sysconfig/ngbauth-submit \
         -B ${XDG_RUNTIME_DIR}  --env KRB5CCNAME="FILE:${XDG_RUNTIME_DIR}/krb5cc"  \
         /cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-analysis/general/pocketcoffea:lxplus-el9-stable


# Create a local virtual environment using the packages defined in the apptainer image
python -m venv --system-site-packages pocket_coffea_env

# Activate the environment
source pocket_coffea_env/bin/activate

# Install in EDITABLE mode
pip install -e .[dev]

cd ../ 
git clone git@github.com:matteomalucchi/PtRegressionJERC.git
cd PtRegressionJERC
pip install -r requirements.txt

# Install the HEPPlotter class
pip install --upgrade  --no-cache-dir git+https://github.com/matteomalucchi/AnalysisConfigs.git
```

After that you should set an alias to activate the PocketCoffea environment because this is called automatically by the `exec.py` script.

On `lxplus`, it can be done by adding the following line to your `~/.bashrc`:

```bash
alias pocket_coffea='apptainer shell --bind /afs -B /cvmfs/cms.cern.ch \
         --bind /tmp  --bind /eos/cms/ -B /etc/sysconfig/ngbauth-submit \
         -B ${XDG_RUNTIME_DIR}  --env KRB5CCNAME="FILE:${XDG_RUNTIME_DIR}/krb5cc"  \
         /cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-analysis/general/pocketcoffea:lxplus-el9-stable'
```

### Other systems installation

If instead you are using a different system, where for example you want to install the environment in micromamba, you can do the following:

```bash
# Clone PocketCoffea
# Clone the fork and checkout the desired branch
git clone https://github.com/PocketCoffea/PocketCoffea.git

# Clone PtRegressionJERC
git git@github.com:matteomalucchi/PtRegressionJERC.git

# Create a local environment and install the packages
cd PtRegressionJERC
micromamba env create -f pocket-coffea-environment.yml
micromamba activate pocket-coffea-env

# Install the HEPPlotter class
pip install --upgrade  --no-cache-dir git+https://github.com/matteomalucchi/AnalysisConfigs.git
```

After that you should set an alias to activate the PocketCoffea environment because this is called automatically by the `exec.py` script.
On your system, it can be done by adding the following line to your `~/.bashrc`:

```bash
alias pocket_coffea='micromamba activate pocket-coffea'
```

### Update HEPPlotter

> [!IMPORTANT]
> To Install the `HEPPlotter` class you can use
>
> ```bash
> pip install --upgrade  --no-cache-dir git+https://github.com/matteomalucchi/AnalysisConfigs.git
> ```
>
> This command should be executed every time you want to pull from the AnalysisConfigs repository and update the `HEPPlotter`.
> If it doesn't update, you should first uninstall it with `pip uninstall configs` and then install it again with the command above.

## Activate the environment

### lxplus environment

To activate the environment, you can use the alias defined above:

```bash
source PocketCoffea/pocket_coffea_env/bin/activate
export PYTHONPATH=$PWD/PocketCoffea:$PYTHONPATH
```

### Other systems environment

To activate the environment, you can use the alias defined above:

```bash
pocket_coffea
```

## JER MC plotting configuration

The script `mc_truth_ptreg_jerc/response_plot/plot_jer_mc.py` is driven by YAML
config files under `mc_truth_ptreg_jerc/response_plot/plot_jer_configs/`. A
config is selected with `-c <path/to/config.yaml>`. Each config has three
groups of settings: top-level flags, bin/variable dictionaries, and per-variable
fit/uncertainty controls.

### Top-level options

| Key | Type | Description |
| --- | --- | --- |
| `mixed_mode` | bool | When true, plot regular and neutrino response variables together using `response_variables_mixed`. Otherwise plot `response_variables` and `response_variables_neutrino` separately. |
| `resolution_vs_pt_gen` | bool | Produce resolution-vs-pT plots using bin centers (`p_T^{ptcl}`). |
| `resolution_vs_pt_reco` | bool | Produce resolution-vs-pT plots using the mean of a mapping variable (e.g. `p_T^{reco}`), and fit the resolution. |
| `histograms_map` | bool | Plot the 1D mapping-variable slices. |
| `histograms_response` | bool | Plot the 1D response-variable slices with optional Gaussian overlay. |
| `gaussian_fit_resolution` | bool | When true, the resolution per bin comes from a Gaussian fit. When false, a confidence-interval (CI) width is used instead. |
| `gaussian_fit_cut_tails` | list/float/false | CI fractions to try when trimming the fit window (e.g. `[0.87]`). `false` disables tail trimming. |
| `gaussian_fit_max_sigma_rel_err` | float | Reject fits where `sigma_err / sigma` is above this value. |
| `min_events_for_fit` | int | Minimum number of events in a 1D slice before attempting a fit / CI. |
| `ci_conf_level` | float | Confidence level used by the CI estimator when Gaussian fits are disabled. |
| `y_lim_resolution` | `[lo, hi]` | Vertical range for the resolution-vs-pT plot. |
| `year_map` | dict | Maps dataset year tags to JERC text-file year tags. |
| `test_*_bins` | list | Reduced bin arrays applied when `--test` is passed on the command line. |

### Bin / mapping / response variable dictionaries

- `bin_variables`, `bin_variables_neutrino` (and `bin_variables_mixed`, built
  automatically as their union if absent): each entry defines an axis used to
  bin response/mapping histograms. Common keys: `bin_edges`, `label`,
  `name_plot`, `resolution_x_variable` (the pT axis used for the resolution
  fit), `txt_name`, `txt_map_to` (for JERC text export).
- `mapping_variables`: variables whose mean per bin is used for the
  resolution-vs-pT plot or to map binning between `mu` and `rho`. Each entry
  needs `bin_vars` plus optional plot parameters (`bin_limits`, `N_bins`,
  `label`, `name_plot`, `legend_name`, `linear_fit`, `txt_name`,
  `rebin_for_plotting`).
- `response_variables`, `response_variables_neutrino`,
  `response_variables_mixed`: variables whose distribution is fitted to extract
  the resolution. See the next section for the keys.

### Per-response-variable options

In addition to the basic keys (`label`, `name_plot`, `N_bins`, `bin_limits`,
`bin_vars`, `legend_name`, `map_x_variable`, `is_reference`, `txt_jet_name`,
`rebin_for_plotting`, `eta_max`), each response variable supports the
following uncertainty / mean controls. Defaults preserve the previous
behaviour.

| Key | Type | Default | Description |
| --- | --- | --- | --- |
| `normalize_by_response_mean` | bool | `false` | When true, plot the relative resolution `sigma / <R>` instead of `sigma`. The uncertainty is propagated as `(sigma/<R>) * sqrt((sigma_err/sigma)^2 + (mean_err/<R>)^2)`. |
| `additional_uncertainty` | float/string | `"0.0"` | Relative uncertainty added in quadrature: `new_err = sqrt(err^2 + (v * sigma)^2)`. Useful to inflate the error band for specific algorithms. |
| `add_min_err` | float/string | `"0.0"` | Minimal uncertainty floor added in quadrature: `new_err = sqrt(err^2 + add_min_err^2)`. For example `"0.0005"` enforces a 0.05% floor. |
| `add_chi2_ndof` | bool | `false` | When true, multiply the uncertainty by `sqrt(chi2 / ndf)` wherever the per-bin Gaussian fit's reduced chi-square exceeds 1. Requires `gaussian_fit_resolution: true`. |
| `response_mean_method` | string | `"auto"` | Selects how the mean of the response is computed per bin. Affects the value used by `normalize_by_response_mean`. Allowed values: `"auto"` (Gaussian fit when available, otherwise binned moments), `"gaussian_fit"` (force the Gaussian fit mean), `"histogram"` (first moment of the binned 1D slice), `"mean_storage"` (mean read from a `hist.storage.Mean()` histogram filled with the un-binned response values; built automatically by `create_ND_histo`). |

The transformations are applied in the order: normalize → `add_min_err` →
`add_chi2_ndof` → `additional_uncertainty`, matching the convention used by
upstream JER tooling.

### Example

```yaml
response_variables:
  MatchedJets_ResponseJEC:
    label: '$p_{\mathrm{T}}^{\mathrm{reco}}$/$p_{\mathrm{T}}^{ptcl}$'
    name_plot: response
    N_bins: *n_bins_response
    bin_limits: *bin_limits_response
    bin_vars: [Pileup_nPU, MatchedJets_eta, MatchedJets_pt]
    legend_name: JEC
    map_x_variable: MatchedJets_JetPtJEC
    is_reference: true
    txt_jet_name: AK4PFPuppi
    rebin_for_plotting: true
    normalize_by_response_mean: true
    additional_uncertainty: "0.2"
    add_min_err: "0.0005"
    add_chi2_ndof: true
    response_mean_method: "mean_storage"
```
