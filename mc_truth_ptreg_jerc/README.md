# MC Truth corrections for pT regression

Repository to compute MC Truth corrections and JER MC for regressed pT jets.

This repository is structured as an analysis configurations for [PocketCoffea](https://github.com/PocketCoffea/PocketCoffea/tree/main) and was originally hosted [here](https://github.com/matteomalucchi/AnalysisConfigs/tree/legacy-jme-config/configs/jme).

## Workflow

### Running the analysis

To run this over the full dataset for a particular year in each $\eta$ and $p_T$ bin, you can use the following command:

```bash
python exec.py --full [-pnet] [-upart] --dir <dir_name> -y <year> [--lxplus]
```

Where `<dir_name>` is the name of the directory where you want to save the results, and `<year>` is the year you want to run the analysis for. The `--lxplus` flag is used to indicate that you are running this on `lxplus` and it will use the `pocket_coffea_env` environment.

Year can be set to:

- 2016_PreVFP
- 2016_PostVFP
- 2017
- 2018
- 2022_preEE
- 2022_postEE
- 2023_preBPix
- 2023_postBPix
- 2024
- 2025

This will save the results in the `dir_name` directory inside the
`output_all.coffea` file.

If running on `lxplus`, there will be an output file for each worker in the `dir_name` directory and you can merge them using:

```bash
cd <dir_name>
pocket-coffea merge-outputs -o output_all.coffea output_job_*.coffea
```

The output file contains 2D histograms for each $\eta$ bin in which the x-axis is the jet $p_T$ response and the y-axis is the jet $p_T$.

#### Running pocket-coffea directly

Instead of `exec.py`, you can also invoke `pocket-coffea` directly. This is useful when you need fine-grained control over the executor, output path, or run options. Set the `YEAR` environment variable beforehand (it is read inside `jer_mc_config.py`):

```bash
export YEAR=2022_postEE
pocket-coffea run \
    --cfg jer_mc_config.py \
    --custom-run-options params/t3_run_options_jer.yaml \
    -o /work/<user>/out_jer_mc_ptreg/out_jer_closure_pnet_upart_2022_postEE \
    --process-separately \
    -e dask@T3_CH_PSI
```

Key flags:

| Flag | Description |
|---|---|
| `--cfg` | PocketCoffea configuration file (`jer_mc_config.py`) |
| `--custom-run-options` | YAML file with executor settings (cores, memory, walltime, chunk size, etc.) |
| `-o` | Output directory where per-dataset coffea/parquet files are written |
| `--process-separately` | Process each dataset independently (enables partial recovery) |
| `-e dask@T3_CH_PSI` | Run on the T3 PSI cluster via Dask; replace with `iterative` for local testing |

The run-options file (`params/t3_run_options_jer.yaml`) controls the Dask workers:

```yaml
scaleout: 30          # number of parallel workers
cores-per-worker: 1
mem-per-worker: "10GB"
disk-per-worker: "1GB"
queue: "standard"
walltime: "08:00:00"
chunksize: 100000
```

Use `params/t3_run_options_jer_15Gb.yaml` if jobs run out of memory on the standard configuration.

The coffea output is written inside the directory specified by `output_chunks_name` in `jer_mc_config.py`. By default this points to the T3 dCache store. The output files are used as input to `plot_jer_mc.py`.

### Computing the MC Truth corrections

After running the full dataset, in order to compute the MC Truth corrections, you can use the following command:

```bash
cd response_plot/
python response.py --full -d <dir_name> --histo 
```

To run on SLURM on tier3:

```bash
cd response_plot/
sbatch -p short --account=t3 --time=00:10:00 --mem 15gb --cpus-per-task=32 --wrap="python response.py --full -d  <dir_name> --histo -n 32"
```

This will:

- Compute the median of the response in each bin in $\eta$ as a function of $p_T$.
- Get the inverse of the median.
- Fit the inverse of the median with a 6th order polynomial.
- Save the results in the configuration file.

It will also:

- Plot the histograms of the response in each bin in $\eta$ and $p_T$ bin.
- Plot the median of the response in each bin in $\eta$ as a function of $p_T$.
- Plot the inverse of the median in each bin in $\eta$ as a function of $p_T$.
- Plot the resolution of the response in each bin in $\eta$ as a function of $p_T$ using 3 different definitions.

### Closure test

To run the closure test of the corrections you can re-run the analysis with some additional flags:

```bash
python exec.py --full -pnet --dir <dir_name> -y <year> --closure --abs-eta-inclusive [--lxplus]
```

This will run the analysis applying the newly derived corrections which have to be specified in the config file.
Once this is done, you can run the other steps of the analysis to obtain the final plots.

To plot all eta bins on the same plot you can use the following command:

```bash
cd response_plot/
python plot_summary_reponse.py -d <dir_name>
```

This is useful to plot the closure test of the MC Truth corrections in a inclusive way.

### JER plotting with `plot_jer_mc.py`

`response_plot/plot_jer_mc.py` is the main plotting script for JER (Jet Energy Resolution) studies from MC truth. It reads the parquet/coffea output produced by the pocket-coffea analysis, builds multi-dimensional histograms of the jet response and auxiliary variables, fits Gaussian distributions to extract the resolution, and produces a full set of publication-quality plots.

#### What the script does

1. **Loads data** — reads all `.coffea` files in the input directory. Each file contains per-chunk arrays of jet variables (response, $p_T$, $\eta$, pile-up, $\rho$, etc.) for each event category and correction type (JEC, Raw, PNet regression, UParT regression, with and without neutrinos).
2. **Builds mapping-variable histograms** — fills N-dimensional histograms for auxiliary variables (`Rho_fixedGridRhoFastjetAll`, reconstructed $p_T$ for each correction) binned in nPU / $\eta$ / gen-$p_T$. A linear fit of $\rho$ vs $\mu$ is performed to convert nPU bin edges into $\rho$ values for the JER text output.
3. **Builds response histograms** — fills per-bin response histograms ($p_T^\text{reco}/p_T^\text{ptcl}$) in a 3-D grid of (nPU, $\eta$, gen-$p_T$). Separate histograms are kept for regular jets and neutrino-inclusive jets.
4. **Fits the response** — fits a Gaussian (with optional tail clipping) to each response histogram to extract the mean and $\sigma$.
5. **Computes JER** — derives the relative resolution $\sigma/\mu$ in each (nPU, $\eta$) slice as a function of gen-$p_T$, and optionally fits the NSC formula: $\sigma/p_T = \sqrt{N^2/p_T^2 + S^2 p_T^d + C^2}$.
6. **Saves results** — writes a coffea file containing all histograms and fit results so they can be reloaded later without re-running the full chain (see `--load` / `--refit`).
7. **Writes JER text files** — exports resolution results in the standard JME text format for each jet correction type and $\eta$ bin.

#### Command-line usage

```bash
cd response_plot/
python plot_jer_mc.py \
    -i <input_dir> \
    -o <output_dir> \
    -c plot_jer_configs/plot_config_jer_mc_coarsePU_rebinResponse.yaml \
    [-w <N_workers>] \
    [--histo] \
    [--novars] \
    [-m <max_parquet>] \
    [-l <precomputed.coffea>] \
    [--refit] \
    [--test]
```

#### Arguments

| Argument | Default | Description |
|---|---|---|
| `-i`/`--input-dir` | required | Directory containing `.coffea` output files from pocket-coffea |
| `-o`/`--output` | `./plot_jer_mc/` | Output directory for all plots, coffea file, and log |
| `-c`/`--config` | required | Path to the YAML configuration file (see below) |
| `-w`/`--workers` | `1` | Number of parallel workers for multiprocessing (speeds up histogram filling and Gaussian fitting) |
| `--histo` | off | Also produce 1D/2D histogram plots of the response and mapping variables |
| `--novars` | off | Use the old save format where variation arrays are not stored inside the coffea files |
| `-m`/`--max-parquet` | `None` (all) | Limit the number of parquet files loaded per dataset (useful for quick tests) |
| `-l`/`--load` | `None` | Path to a precomputed coffea file produced by a previous run; skips data loading and histogram filling |
| `--refit` | off | When used with `--load`, re-merge histogram bins according to the current YAML config (coarser binning) and recompute all Gaussian fits before plotting |
| `--test` | off | Run in test mode with reduced bin arrays (overrides the `test` flag in the YAML config) |

#### Typical workflows

**First run from scratch:**

```bash
python plot_jer_mc.py -i /work/<user>/out_jer_mc_ptreg/out_jer_closure_pnet_upart_2022_postEE \
    -o plots_2022postEE \
    -c plot_jer_configs/plot_config_jer_mc_coarsePU_rebinResponse.yaml \
    -w 16 --histo
```

**Reload and replot with coarser binning (no re-reading of parquet):**

```bash
python plot_jer_mc.py \
    -l plots_2022postEE/plotted_data.coffea \
    -o plots_2022postEE_rebin \
    -c plot_jer_configs/plot_config_jer_mc_coarsePU_rebinResponse.yaml \
    --refit -w 16
```

**Quick test on a subset of files:**

```bash
python plot_jer_mc.py -i /work/<user>/out_jer_... -o test_out \
    -c plot_jer_configs/plot_config_jer_mc_coarsePU_rebinResponse.yaml \
    -m 5
```

(Use `test: true` in the YAML config to also reduce the bin granularity.)

---

### Plotting JERC txt files with `plot_jerc_txt.py`

`response_plot/plot_jerc_txt.py` reads one or more JME-format resolution (or L2Relative correction) txt files and produces publication-quality plots of the resolution (or jet response) as a function of jet $\eta$, evaluated at a set of user-chosen $p_T$ values.

#### What the script does

1. **Parses the JERC txt format** — reads the file header (bin variables, x variable, formula string) and each row of bin edges and formula parameters. Rows with non-finite parameters are logged and skipped.
2. **Evaluates the formula** — for every combination of (JetPt, $\eta$ bin, $\rho$ bin), the ROOT `TFormula` string is evaluated numerically via numpy. Unphysical combinations (where $E_\text{jet} = p_T \cosh|\eta| \geq E_\text{beam}/2$) are silently dropped.
3. **Produces resolution-vs-$\eta$ plots** — one plot per input file, with:
   - **Color + marker** per $p_T$ evaluation value.
   - **Alpha shading** per $\rho$ bin (lighter = lower $\rho$, darker = higher $\rho$). When no $\rho$ binning is present a single opaque marker per $p_T$ is used.
   - **CMS detector-region lines and labels** (Barrel, EC1, EC2, HF) overlaid on the $\eta$ axis.
   - A full-range plot (all $\eta$) and a restricted-range plot ($|\eta| < $ `--eta-max`, default 2.5) are saved for each file.
4. **Handles L2Relative correction files** — when the txt header declares `Correction L2Relative`, the script inverts the correction ($1/c$) and plots it as the simulated jet response with a horizontal reference line at 1.
5. **Produces pairwise ratio plots** — for every pair of input files, plots the point-wise ratio (numerator / denominator) as a function of $\eta$ with a horizontal line at 1. Ratio plots are also saved in the restricted-$\eta$ range. Use `--no-ratio` to skip this step.

#### Command-line usage

```bash
cd response_plot/

# single file
python plot_jerc_txt.py Run3Summer22_V1_NSC_MC_PtResolution_AK4PFPuppi.txt \
    --pt-values 80 150 300 600 -o plots/

# all txt files in a directory (shell expands the glob)
python plot_jerc_txt.py /path/to/txt/*.txt --pt-values 80 150 300 600 -o plots/
```

#### Arguments

| Argument | Default | Description |
|---|---|---|
| `input` | required | One or more JERC resolution/correction txt files; shell globs (`*.txt`) are accepted |
| `--pt-values` | `50 100 300 600 1000 3000` | JetPt [GeV] values at which to evaluate the resolution formula |
| `-o`/`--output` | `.` | Output directory for all plots |
| `--formats` | `png pdf svg` | Output file formats |
| `--no-ratio` | off | Skip pairwise ratio plots |
| `--eta-max` | `2.5` | $|\eta|$ upper limit for the restricted-range plots |

#### Output files

For each input file `<stem>.txt` the script writes:

- `<output>/<stem>.<fmt>` — full-$\eta$ resolution (or response) plot.
- `<output>/<stem>_abseta_lt_<eta_max>.<fmt>` — restricted-$\eta$ version.

For each pair of files `<stem_a>` and `<stem_b>`:

- `<output>/ratio_<stem_a>__over__<stem_b>.<fmt>` — full-$\eta$ ratio plot.
- `<output>/ratio_<stem_a>__over__<stem_b>_abseta_lt_<eta_max>.<fmt>` — restricted-$\eta$ ratio.

---

### Plot configuration file (`plot_jer_configs/`)

The YAML configuration file controls every aspect of the binning, variables, response histograms, and plot style. Several pre-made configurations exist in `response_plot/plot_jer_configs/`.

#### Boolean flags

| Key | Default | Description |
|---|---|---|
| `test` | `false` | Replace all bin arrays with the reduced `test_*` arrays for fast iteration |
| `mixed_mode` | `true` | Process regular and neutrino-inclusive response variables together in a single `mixed` mode; when `false` they are processed as two separate `regular` and `neutrino` modes |
| `resolution_vs_pt_gen` | `false` | Plot JER as a function of gen-$p_T$ |
| `resolution_vs_pt_reco` | `false` | Plot JER as a function of mean reco-$p_T$ (mapped x-axis) |
| `histograms_map` | `false` | Save/plot the mapping-variable histograms |
| `histograms_response` | `false` | Save/plot the 1D response histograms per bin |
| `gaussian_fit_resolution` | `true` | Derive resolution from a Gaussian fit to the response histogram |
| `gaussian_fit_cut_tails` | `true` | Clip the response histogram to a symmetric window around the peak before fitting |

#### Bin arrays

```yaml
pu_bins: [0, 20, 30, 40, 50, 60, 100]   # nPU bin edges
jet_eta_bins: [-5.191, ..., 5.191]       # η bin edges (JME standard)
jet_pt_bins: [8, 10, ..., 5000]          # gen-pT bin edges [GeV]
```

Reduced arrays used when `test: true`:

```yaml
test_pu_bins: [0, 50, 100]
test_jet_eta_bins: [-5.191, 0.0, 5.191]
test_jet_pt_bins: [8, 50, 120, ...]
```

#### Response histogram settings

```yaml
n_bins_response: 150          # number of bins in the response (pT_reco/pT_gen) histogram
bin_limits_response: [0.5, 1.5]  # [min, max] of the response axis
y_lim_resolution: [0, 0.35]   # y-axis range for resolution plots
```

#### `year_map`

Maps the year tag extracted from dataset names to the campaign string used in JME text-file names:

```yaml
year_map:
  "2022_preEE": Summer22
  "2022_postEE": Summer22EE
  "2023_preBPix": Summer23
  ...
```

#### `bin_variables` and `bin_variables_neutrino`

Define the variables used as binning axes for the response histograms. Each entry is a column name from the coffea arrays:

```yaml
bin_variables:
  Pileup_nPU:
    event_var: true          # true → one value per event; false → one value per jet
    bin_edges: *pu_bins      # YAML anchor referencing the bin array above
    label: '$\mu$'           # axis label for plots
    name_plot: nPU           # short name used to build plot filenames
    resolution_x_variable: false  # true → this is the x-axis of the JER vs pT plot
    txt_map_to: Rho_fixedGridRhoFastjetAll  # variable to map this axis to in txt output
  MatchedJets_eta:
    event_var: false
    bin_edges: *jet_eta_bins
    label: '$\eta$'
    name_plot: jet_eta
    resolution_x_variable: false
    txt_name: JetEta         # name used in the JME text-file header
  MatchedJets_pt:
    event_var: false
    bin_edges: *jet_pt_bins
    label: '$p_{\mathrm{T}}^{ptcl}$ [GeV]'
    name_plot: jet_gen_pt
    resolution_x_variable: true   # this axis becomes the x-axis of the JER plot
    txt_name: JetPt
```

`bin_variables_neutrino` mirrors `bin_variables` but for the neutrino-inclusive jet collection (`MatchedJetsNeutrino_*`).

#### `mapping_variables`

Variables that are profiled and fit (e.g. $\rho$ vs $\mu$, reco-$p_T$ vs gen-$p_T$). Used to map bin centres to mean values for the x-axis of plots and for the JME text output:

```yaml
mapping_variables:
  Rho_fixedGridRhoFastjetAll:
    event_var: true
    label: '$\rho$ [GeV/area]'
    name_plot: rho
    N_bins: 1000
    bin_vars: [Pileup_nPU]    # bins in which the profile is computed
    linear_fit: true          # fit a linear model ρ = m·μ + b per η/pT bin
    txt_name: Rho
    rebin_for_plotting: true
  MatchedJets_JetPtJEC:
    label: '$p_{\mathrm{T}}^{\mathrm{reco}}$ [GeV]'
    name_plot: jet_reco_pt
    N_bins: 1000
    bin_limits: [0, 6000]
    bin_vars: [Pileup_nPU, MatchedJets_eta, MatchedJets_pt]
    legend_name: JEC
    rebin_for_plotting: true
```

Add or remove entries to include additional correction types or auxiliary variables.

#### `response_variables`, `response_variables_neutrino`, `response_variables_mixed`

Define which response columns to read and how to plot them. Each entry corresponds to one correction type:

```yaml
response_variables:
  MatchedJets_ResponseJEC:
    label: '$p_{\mathrm{T}}^{\mathrm{reco}}$/$p_{\mathrm{T}}^{ptcl}$'
    name_plot: response_jec
    N_bins: *n_bins_response
    bin_limits: *bin_limits_response
    bin_vars: [Pileup_nPU, MatchedJets_eta, MatchedJets_pt]
    legend_name: JEC          # label in the legend
    map_x_variable: MatchedJets_JetPtJEC  # mapping variable for the reco-pT x-axis
    is_reference: true        # marks this as the reference curve in overlay plots
    txt_jet_name: AK4PFPuppi  # jet algorithm string in the JME text-file name
    rebin_for_plotting: true  # rebin the stored fine-binned histogram before plotting
```

`response_variables_mixed` combines regular and neutrino-inclusive entries and is used when `mixed_mode: true`.

#### `plot_settings`

Controls the color and marker style for each correction type in overlay plots. Keys are matched against the suffix of the response variable name:

```yaml
plot_settings:
  JEC:              { color: darkorange, fmt: "o" }
  Raw:              { color: green,      fmt: "s" }
  PNetReg:          { color: darkred,    fmt: "<" }
  UparTReg:         { color: tomato,     fmt: "^" }
  PNetRegNeutrino:  { color: darkblue,   fmt: ">" }
  UparTRegNeutrino: { color: cornflowerblue, fmt: "v" }
```
