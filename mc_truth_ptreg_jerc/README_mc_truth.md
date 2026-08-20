# MC truth corrections from dumped columns

This document describes the **column-based MC truth workflow**: the derivation of
the L2Relative (MC truth) corrections for the regressed jet pT starting from
flat columns dumped into parquet files, in the same spirit as the JER MC study
(`jer_mc_config.py` + `response_plot/plot_jer_mc.py`).

It replaces the histogram-based chain (`mc_truth_config.py` /
`cartesian_config.py` + `exec.py` + `response_plot/response.py`), which is kept
in the repository for reference but is no longer needed. See
[`README.md`](README.md) for the legacy instructions and the description of
`plot_jer_mc.py` / `plot_jerc_txt.py`, and the [top-level
README](../README.md) for the installation.

## What changed

| Old workflow | New workflow |
|---|---|
| One job per eta region (`neg1`, `neg2`, ..., `pos4`) and per flavour | **One job** covering the full eta range and all flavours |
| One job for ParticleNet, another one for UParT | **PNet and UParT together**, with and without neutrinos |
| Output: one 2D histogram per eta/pT category inside `output_all.coffea` | Output: **flat columns dumped into parquet files** |
| Binning frozen at analysis time (`CartesianSelection`, `MultiCut`) | Binning chosen **offline**, in a YAML config, and changed without re-running |
| ~15 environment variables exported in bash (`YEAR`, `SIGN`, `PNET`, `UPART`, `CLOSURE`, `CARTESIAN`, `FLAVSPLIT`, ...) | A single **YAML run card**, written automatically from the command line |
| Jobs launched inside `tmux` sessions by `exec.py` | A plain foreground command (`run_mc_truth.py`); use `nohup`, `screen` or `sbatch` if you want it detached |
| Flavour splitting done at analysis level (one job per flavour) | Flavour splitting done **at plotting level**, from the flavour column |

## New files

| File | Purpose |
|---|---|
| `mc_truth_columns_config.py` | PocketCoffea configuration: selects the dataset, applies the corrections and declares the columns dumped to parquet |
| `mc_truth_workflow.py` | `MCTruthProcessor`: jet matching and response computation for JEC, Raw, PNet and UParT (with and without neutrinos). No environment variables, no per-eta collections |
| `mc_truth_corrections.py` | Environment-free parser of the JME `L2Relative` txt files |
| `params/mc_truth_runcard.yaml` | Default run card: year, dataset, physics options |
| `params/t3_run_options_mc_truth.yaml` | Executor options for the PSI T3 |
| `params/lxplus_run_options_mc_truth.yaml` | Executor options for lxplus/HTCondor |
| `run_mc_truth.py` | Command-line launcher (replaces `exec.py`) |
| `response_plot/response_mc_truth.py` | Derivation of the corrections and all the plots (replaces `response_plot/response.py`) |
| `response_plot/column_utils.py` | Shared machinery: column loading, N-D histograms, per-bin estimators, fits |
| `response_plot/write_l2rel_mc_truth.py` | Writer of the `L2Relative` txt files (replaces `response_plot/write_l2rel.py`) |
| `response_plot/mc_truth_configs/mc_truth_config_default.yaml` | Default plotting/derivation configuration, fully commented |

---

## 1. Run the analysis

```bash
python run_mc_truth.py -y 2023_preBPix -o /work/$USER/out_mc_truth_2023_preBPix
```

That is the whole command: no `export` before it, no tmux session around it.
It writes

* `<output>/mc_truth_runcard.yaml` — the exact configuration used for this run,
* `<output>/t3_run_options_mc_truth.yaml` — the executor options actually used,
* `<output>/output_all.coffea` (or one file per job) — the bookkeeping file,
* the per-chunk **parquet files** in `columns_output_dir` (see below).

### Options

| Option | Default | Description |
|---|---|---|
| `-y`, `--year` | from the run card | Campaign: `2016_PreVFP`, `2016_PostVFP`, `2017`, `2018`, `2022_preEE`, `2022_postEE`, `2023_preBPix`, `2023_postBPix`, `2024`, `2025` |
| `-o`, `--output` | required | Output directory |
| `-c`, `--runcard` | `params/mc_truth_runcard.yaml` | Base run card; the options below are applied on top of it |
| `--columns-dir` | from the run card | Destination of the parquet files; `{year}` is substituted |
| `-e`, `--executor` | `dask@T3_CH_PSI` | Any PocketCoffea executor (`iterative`, `dask@T3_CH_PSI`, `condor@lxplus`, ...) |
| `-r`, `--run-options` | `params/t3_run_options_mc_truth.yaml` | Executor options (workers, memory, walltime, chunk size) |
| `--closure` | off | Apply the freshly derived MC truth corrections on top of the regressions (closure test) |
| `--no-pnet` / `--no-upart` | off | Skip one of the two taggers (by default **both** run) |
| `--no-neutrinos` | off | Skip the neutrino-inclusive gen jet collection |
| `--no-reapply-jec` | off | Use the jet pT calibrated by PocketCoffea instead of re-applying the L2Relative txt file |
| `--test` | off | `pocket-coffea run --test`: a couple of chunks, locally |
| `--process-separately` | off | Process each dataset independently |
| `--log` | off | Tee the output to `<output>/run_mc_truth.log` |
| `--overwrite` | off | Run even if `<output>/output_all.coffea` exists |
| `-n`, `--dry-run` | off | Only write the run card and print the command |

### Running on the T3 batch system

`run_mc_truth.py` runs in the foreground, so submit it the way you prefer:

```bash
sbatch -p standard --account=t3 --time=08:00:00 --mem 8gb \
    --wrap="python run_mc_truth.py -y 2023_preBPix -o /work/$USER/out_mc_truth_2023_preBPix"
```

### Running on lxplus

```bash
python run_mc_truth.py -y 2023_preBPix -o /afs/cern.ch/work/$USER/out_mc_truth_2023_preBPix \
    -e condor@lxplus -r params/lxplus_run_options_mc_truth.yaml
```

HTCondor workers re-import the configuration file, so `run_mc_truth.py` appends
one `custom-setup-commands` line to the copied run-options file that exports
`MC_TRUTH_RUNCARD` with the path of the generated run card. The output directory
must therefore be on a shared file system (AFS/EOS) when running on condor.

### Calling `pocket-coffea` directly

```bash
pocket-coffea run --cfg mc_truth_columns_config.py \
    -o /work/$USER/out_mc_truth_2023_preBPix \
    -e dask@T3_CH_PSI --custom-run-options params/t3_run_options_mc_truth.yaml
```

Without `MC_TRUTH_RUNCARD` the configuration reads
`params/mc_truth_runcard.yaml`, so edit that file (or copy it and export
`MC_TRUTH_RUNCARD=/path/to/your_runcard.yaml`) to change the year or any other
option.

### The run card

`params/mc_truth_runcard.yaml` is fully commented; the keys are:

| Key | Description |
|---|---|
| `year` | Campaign to process |
| `dataset_jsons`, `samples` | Dataset definition files and the sample name per year |
| `columns_output_dir` | Where the parquet chunks go (`{year}` is substituted) |
| `pnet`, `upart`, `neutrinos` | Which responses are computed |
| `deltaR_matching` | Cone used to match reco jets to gen jets |
| `genjet_pt_cut` | Minimum gen jet pT (0: the pT binning is applied offline) |
| `n_jets` | Number of leading raw-pT jets kept per event |
| `com_energy` | Centre-of-mass energy, used to drop unphysical jets |
| `reapply_jec` | Re-apply the L2Relative txt file on the raw pT to define the JEC response |
| `closure` | Apply the derived corrections on top of the regressions |
| `regression_value_when_invalid` | `raw`, `original` or a float, used when a regression factor is not positive |
| `invalid_regression_eta` | Above this \|eta\| the regression falls back to the raw jet |
| `pv_presel_distance` | Maximum distance between the reconstructed and the generated primary vertex |

### Dumped columns

For `MatchedJets` (gen jets) and `MatchedJetsNeutrino` (gen jets including the
matched neutrinos):

`pt` (gen), `eta` (gen), `RecoEta`, `RecoPhi`, `partonFlavour`, `hadronFlavour`,
`ResponseJEC`, `JetPtJEC`, `ResponseRaw`, `JetPtRaw`, `ResponsePNetReg`,
`JetPtPNetReg`, `ResponseUparTReg`, `JetPtUparTReg` (and the corresponding
`...Neutrino` versions), plus `Rho_fixedGridRhoFastjetAll` and `Pileup_nPU`.

`partonFlavour` and `hadronFlavour` are what makes the flavour splitting
possible at plotting time.

---

## 2. Derive the corrections and make the plots

```bash
cd response_plot/
python response_mc_truth.py -i /work/$USER/out_mc_truth_2023_preBPix \
    -o /work/$USER/out_mc_truth_2023_preBPix/mc_truth_plots -w 16
```

The script:

1. reads the parquet columns and flattens them to one entry per jet;
2. assigns each jet to a flavour group and fills one N-dimensional histogram per
   response variable, with axes `(flavour, eta^reco, pT^ptcl, response)`;
3. computes per bin the **median** of the response, its **inverse** (the
   correction), the **resolution** and the **mean reconstructed pT**;
4. fits the inverse median with a polynomial in `log10(pT^reco)`, increasing the
   order until the p-value exceeds `fit_p_value_threshold`;
5. writes the `L2Relative` txt files, inclusively and per flavour;
6. produces the plots.

Because the flavour is an axis of the histogram and every jet belongs to exactly
one group, the inclusive result is the sum over that axis — i.e. exactly the
sample you would get without any splitting.

### Options

| Option | Default | Description |
|---|---|---|
| `-i`, `--input-dir` | — | Directory with the `.coffea` files written by `run_mc_truth.py` |
| `-o`, `--output` | `./mc_truth_plots` | Output directory |
| `-c`, `--config` | built-in default | YAML config merged on top of `mc_truth_configs/mc_truth_config_default.yaml` |
| `-w`, `--workers` | `1` | Processes used to draw the plots |
| `-m`, `--max-parquet` | all | Maximum number of parquet files read per dataset (quick tests) |
| `-l`, `--load` | — | Replot from a previously saved `mc_truth_results_<category>.coffea` |
| `--histo` | off | Also plot the 1D response distribution of every bin |
| `--max-histo-plots` | `2000` | Safety limit on the number of response histogram plots |
| `--no-plots` | off | Only derive the corrections (fits, JSON, txt) |
| `--no-save-histograms` | off | Do not store the response histograms in the results file (smaller output, but `--load` cannot redraw them) |
| `--no-flavour-split` | off | Process the inclusive flavour only |
| `--novars` | off | Old save format, without the variations in the coffea files |
| `--test` | off | Use the reduced eta/pT test binning |

On SLURM:

```bash
cd response_plot/
sbatch -p short --account=t3 --time=01:00:00 --mem 30gb --cpus-per-task=32 \
    --wrap="python response_mc_truth.py -i <input_dir> -o <output_dir> -w 32 --histo"
```

### Outputs

```
<output>/
├── mc_truth_config_default.yaml          # copy of the configuration used
├── response_mc_truth.log
├── fit_results_inverse_median_baseline.json   # every fit: parameters, chi2, p-value, points
├── mc_truth_results_baseline.coffea      # everything, reloadable with --load
├── l2relative_txt/
│   ├── Summer23Run3_V3_MC_L2Relative_AK4PFPNet.txt
│   ├── Summer23Run3_V3_MC_L2Relative_AK4PFPNet_bFlav.txt        # if requested
│   ├── Summer23Run3_V3_MC_L2Relative_AK4PFPNetPlusNeutrino.txt
│   ├── Summer23Run3_V3_MC_L2Relative_AK4PFUparT.txt
│   └── Summer23Run3_V3_MC_L2Relative_AK4PFUparTPlusNeutrino.txt
└── plots_baseline/
    ├── median/            # median response vs gen pT, one plot per (flavour, eta bin)
    ├── inv_median/        # 1/median vs mean reco pT, with the polynomial fit
    ├── resolution/        # resolution vs gen pT, with the ratio to JEC
    ├── flavour_median/    # flavours overlaid, one plot per (correction, eta bin)
    ├── flavour_resolution/
    └── response_histograms/   # only with --histo
```

The txt files can be copied into `params/` and plotted/compared with
`plot_jerc_txt.py` (see [`README.md`](README.md)).

### Replotting without re-reading the parquet files

```bash
python response_mc_truth.py -l <output>/mc_truth_results_baseline.coffea -o <output>_replot
```

Everything that depends only on the already-computed numbers (all the plots) is
redone. Changing the **binning** requires re-running from the columns, which is
cheap compared to re-running the analysis.

---

## 3. Configuration of the derivation (`mc_truth_configs/`)

`mc_truth_configs/mc_truth_config_default.yaml` is the built-in default and is
annotated option by option. A config passed with `-c` is merged on top of it:
scalars and lists are replaced, dictionaries are merged key by key, so
overriding a single option of a single variable is enough:

```yaml
# my_config.yaml — coarser eta binning and the confidence-interval resolution
jet_eta_bins: &jet_eta_bins [-5.191, -3.0, -2.5, -1.3, 0.0, 1.3, 2.5, 3.0, 5.191]
resolution_estimator: confidence
bin_variables:
  MatchedJets_RecoEta:
    bin_edges: *jet_eta_bins
  MatchedJetsNeutrino_RecoEta:
    bin_edges: *jet_eta_bins
```

### Main options

| Key | Default | Description |
|---|---|---|
| `test` | `false` | Use the reduced `test_jet_eta_bins` / `test_jet_pt_bins` |
| `min_events_per_bin` | `20` | Bins with fewer jets give NaN |
| `median_min` | `0.1` | Bins whose median is below this are considered broken and dropped |
| `resolution_estimator` | `quantile` | `quantile` ((p84 − p16)/2) or `confidence` (smallest interval holding `ci_conf_level`) |
| `ci_conf_level` | `0.87` | Confidence level of the `confidence` estimator |
| `normalize_resolution_by_median` | `true` | Plot sigma/median instead of the bare width |
| `rebin_for_plotting` | `true` | Rebin the response histogram before drawing the 1D slices |
| `num_params` | `9` | Columns per row of the txt file: parameters + 2. `9` means at most a pol6 |
| `fit_p_value_threshold` | `0.05` | The polynomial order grows until the p-value exceeds this |
| `fit_pt_min` | `1.0` | Only bins with a mean reco pT above this enter the fit |
| `write_txt` | `true` | Write the `L2Relative` txt files |
| `version` | `V3` | Version tag used in the txt file names |
| `campaign_map` | — | Year tag → campaign name used in the txt file names |
| `n_bins_response` | `500` | Bins of the response axis (see the note on memory in the file) |
| `bin_limits_response` | `[0, 3]` | Range of the response axis |
| `jet_eta_bins`, `jet_pt_bins` | JME standard | Binning of the derivation |
| `y_lim_median`, `y_lim_inverse_median`, `y_lim_resolution` | — | y ranges of the summary plots |

The median and the quantiles are obtained by interpolating the cumulative
distribution inside the bin, so a moderate `n_bins_response` is enough: with 500
bins over `[0, 3]` the median agrees with the un-binned one at the 1e-5 level.

### `bin_variables`, `mapping_variables`, `response_variables`

Same structure as `plot_jer_configs/plot_config_jer_mc_default.yaml`:

* **`bin_variables`** are the axes of the histograms. `kind: flavour` marks the
  categorical flavour axis (its categories are filled in automatically),
  `x_variable: true` marks the axis used as x-axis of the summary plots (the gen
  pT), and the remaining axis is the eta axis. `txt_name` is the name used in the
  txt file header.
* **`mapping_variables`** are profiled per bin; their mean is used as x-axis of
  the inverse median plot and, most importantly, as the **x variable of the
  fit**: the correction is a function of the *reconstructed* pT. `mean_only:
  true` stores only the per-bin mean and skips the count histogram, which saves a
  lot of memory.
* **`response_variables`** are the responses themselves. `derive_correction:
  true` requests the polynomial fit and the txt file, whose name uses
  `txt_jet_name`; `is_reference: true` marks the denominator of the ratio panel;
  `eta_max` skips the bins outside the training region of a tagger (2.5 for
  UParT).

---

## 4. Flavour splitting

The flavour splitting is now a plotting-time operation: no dedicated job, no
`FLAVSPLIT`/`FLAV` environment variables.

```yaml
flavour_split:
  enabled: true
  field: partonFlavour        # or hadronFlavour
  use_abs: true               # compare |pdgId|
  groups:
    uds: [1, 2, 3]
    c: [4]
    b: [5]
    g: [21]
  other_name: other           # catch-all: everything that matches no group
  plot_flavours: [inclusive, b, c, uds, g]
  derive_correction_flavours: [inclusive]
  labels:  {inclusive: all flavours, b: b, c: c, uds: uds, g: gluon}
  colors:  {inclusive: black, b: darkblue, c: darkgreen, uds: darkorange, g: darkred}
  markers: {inclusive: "o", b: "s", c: "^", uds: "v", g: "<"}
```

* `plot_flavours` — flavours shown in the plots. `inclusive` is the sum over the
  whole flavour axis, so it is exactly the un-split sample.
* `derive_correction_flavours` — flavours for which the inverse median is fitted
  and a txt file is written. Flavour-specific files get a `_<flav>Flav` suffix,
  e.g. `Summer23Run3_V3_MC_L2Relative_AK4PFPNet_bFlav.txt`.
* Every jet belongs to exactly one group; jets whose flavour matches no group
  (e.g. `partonFlavour == 0`) end up in `other_name`, so nothing is lost and the
  sum over the axis is the inclusive sample.

Two families of plots come out of it:

* `plots_<category>/median/`, `inv_median/`, `resolution/`: for a **given
  flavour**, all the corrections overlaid (one plot per eta bin);
* `plots_<category>/flavour_median/`, `flavour_resolution/`: for a **given
  correction**, all the flavours overlaid, with the ratio to the inclusive
  sample (one plot per eta bin).

Use `--no-flavour-split` to switch the whole machinery off (the flavour axes are
then removed from the histograms, which also makes the job lighter).

---

## 5. Closure test

Derive the corrections, copy the txt files where the configuration expects them,
and re-run with `--closure`:

```bash
# 1. derive
python run_mc_truth.py -y 2023_preBPix -o /work/$USER/out_mc_truth_2023_preBPix
cd response_plot/
python response_mc_truth.py -i /work/$USER/out_mc_truth_2023_preBPix \
    -o /work/$USER/out_mc_truth_2023_preBPix/mc_truth_plots -w 16
cd ..

# 2. install the new txt files (the names are listed in params/l2relative_txt.py)
cp /work/$USER/out_mc_truth_2023_preBPix/mc_truth_plots/l2relative_txt/*.txt params/

# 3. re-run applying them, and derive again: the medians must now sit at 1
python run_mc_truth.py -y 2023_preBPix -o /work/$USER/out_mc_truth_2023_preBPix_closure --closure
cd response_plot/
python response_mc_truth.py -i /work/$USER/out_mc_truth_2023_preBPix_closure \
    -o /work/$USER/out_mc_truth_2023_preBPix_closure/mc_truth_plots -w 16
```

`params/l2relative_txt.py` maps each year to the txt file names that are picked
up, both for the JEC re-application and for the closure test.

---

## 6. Relation to `plot_jer_mc.py`

The two analyses now share the same input format (columns in parquet), the same
YAML-driven binning and the same `HEPPlotter` output. The generic pieces
(configuration merging, column flattening, N-D histogram filling, automatic
rebinning, per-bin estimators, fits) live in `response_plot/column_utils.py`,
which `response_mc_truth.py` is built on; `plot_jer_mc.py` still carries its own
copies of the equivalent helpers and can be migrated to `column_utils.py` later.

In practice:

* **JER MC** (`jer_mc_config.py` + `plot_jer_mc.py`): resolution from a Gaussian
  fit, binned in nPU/eta/pT, NSC fit, `PtResolution` txt files.
* **MC truth** (`mc_truth_columns_config.py` + `response_mc_truth.py`): median
  response, binned in flavour/eta/pT, polynomial fit of the inverse median,
  `L2Relative` txt files.
