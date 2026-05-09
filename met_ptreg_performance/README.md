# MET Studies for regressed pT jets

Repository to compute the response and resolution of MET as well as the effect of the Type-1 correction derived for PNet regressed pT jets.

This repository is structured as an analysis configuration for [PocketCoffea](https://github.com/PocketCoffea/PocketCoffea/tree/main) and was originally hosted [at this link](https://github.com/matteomalucchi/AnalysisConfigs/tree/legacy-met-config/configs/MET_studies).

## What the program does

### Physics overview

The analysis uses $Z\to\mu\mu$ events as a clean probe to measure MET performance. Since the Z boson transverse momentum $\vec{q}_T$ is precisely known from the two muons, it can be balanced against the hadronic recoil $\vec{u}$ to quantify how well MET is reconstructed. A perfect detector would give $\vec{u} = -\vec{q}_T$.

### Event processing ([`workflow.py`](workflow.py))

The `METProcessor` class performs the following steps:

1. **Jet preparation** — Several jet collections are built from NanoAOD `Jet` and `CorrT1METJet` branches:
   - Standard JEC-corrected jets (`JetMuonSubtr`)
   - PNet-regressed pT jets (`JetPNetMuonSubtr`)
   - PNet-regressed pT + neutrino correction jets (`JetPNetPlusNeutrinoMuonSubtr`)
   - Low-pT jets from the `CorrT1METJet` collection (`JetLowPtMuonSubtr`), used to recover soft hadronic activity below the standard jet pT threshold
   - All jet collections have muon four-momenta subtracted (muon subtraction factor and optionally `muonSubtrDeltaEta`/`muonSubtrDeltaPhi` angle corrections) to avoid double-counting with the well-measured muons

2. **Type-1 MET correction** — Raw PuppiMET is propagated to several corrected MET branches by replacing raw jet $p_T$ with JEC-corrected or PNet-regressed $p_T$ (`RawPuppiMET-Type1<suffix>`). The correction is applied as:

$$\vec{p}_T^{\,\text{MET,corr}} = \vec{p}_T^{\,\text{RawPuppiMET}} - \sum_\text{jets} \vec{p}_T^{\,\text{raw}} + \sum_\text{jets} \vec{p}_T^{\,\text{corr}}$$

   Multiple options (`option_1` through `option_6`) control how jets without a valid regression output are handled (masked out, replaced by JEC jets, etc.).

1. **Lepton and dilepton selection** — Muons and electrons are selected; the two leading muons form the Z candidate (`ll`), whose transverse momentum defines $\vec{q}_T$.

2. **Hadronic recoil and projections** — For each MET branch the hadronic recoil is computed as:

$$\vec{u} = -\left(\vec{p}_T^{\,\text{MET}} + \vec{q}_T\right)$$

   and projected onto two components relative to the Z direction:

- $u_\parallel$ (`u_paral`): component of $\vec{u}$ parallel to $\vec{q}_T$, shifted by $|\vec{q}_T|$ so that a perfect response gives $u_\parallel = 0$
- $u_\perp$ (`u_perp`): component of $\vec{u}$ perpendicular to $\vec{q}_T$, which should be zero on average
- $R$ (response): scalar response defined as $R = \vec{u}\cdot(-\hat{q}_T)\,/\,|\vec{q}_T|$, equal to 1 for perfect MET

### Summary plots ([`plot_MET.py`](plot_MET.py))

The plotting script reads the coffea output files and produces three categories of plots, all comparing the different MET types (standard Type-1, PNet-regressed Type-1, etc.) simultaneously:

1. **Inclusive MET histograms** — pT distributions of the various MET collections, including relative difference plots between standard PuppiMET and the Type-1 corrected variants.

2. **Response and resolution curves** — For each MET type the following quantities are computed in bins of $|\vec{q}_T|$ and number of primary vertices ($n_\text{PV}$):
   - Mean response $\langle R \rangle$ (ideal value: 1)
   - Mean $u_\parallel$ and $u_\perp$ (ideal value: 0)
   - Quantile resolution: $(Q_{84} - Q_{16})/2$ of $u_\parallel$ and $u_\perp$
   - Standard deviation resolution of $u_\parallel$ and $u_\perp$
   - Scaled variants $u_\parallel / \langle R \rangle$ and $u_\perp / \langle R \rangle$ (divided by $\langle R \rangle$ in each bin) to disentangle response and resolution effects

3. **1D and 2D response histograms** — Per-bin distributions of $R$, $u_\parallel$, $u_\perp$ (and their scaled versions) as 1D histograms, and 2D histograms of response vs binning variable with profile plots of the mean and resolutions. These are produced only when the `--histo` flag is passed.

## Workflow

### Running the analysis

To run the analysis on Tier3, use the following command:

```bash
pocket-coffea run --cfg MET_studies_config.py --custom-run-options params/t3_run_options_big.yaml -o <output-dir> -e dask@T3_CH_PSI
```

To produce the response plots, use:

```bash
python plot_MET.py -i <input-dir> -w 8 --histo -o <output-plot-dir>
```

### Lxplus

The easiest way to run on lxplus is via the provided script, which handles venv activation and `PYTHONPATH` setup automatically. Enter the singularity first, then run it from inside:

```bash
# activate singularity
pocket_coffea

# run the script (output dir is required)
cd PtRegressionJERC/met_ptreg_performance
./run_lxplus.sh <output-dir>
```

Alternatively, the same steps can be run manually:

```bash
# activate singularity
pocket_coffea

# activate venv
cd PocketCoffea
source pocket_coffea_env/bin/activate
# Set the PYTHONPATH to make sure the editable PocketCoffea installation is picked up
export PYTHONPATH=`pwd`

# run pocket coffea
cd ../PtRegressionJERC/met_ptreg_performance 
pocket-coffea run --cfg MET_studies_config_lxplus.py --custom-run-options ./params/lxplus_run_options_big.yaml -o <output-dir> -e condor@lxplus
```

> [!IMPORTANT]
> Before running, update the `output_chunks_name` path in `MET_studies_config_lxplus.py` to point to your own EOS area:
>
> ```python
> output_chunks_name = "root://eosuser.cern.ch//eos/user/<initial>/<username>/..."
> ```

### Adding new datasets

New datasets are defined in [`datasets/datasets_definitions.json`](datasets/datasets_definitions.json). Each entry follows the PocketCoffea format:

```json
{
    "MySample": {
        "sample": "MySample",
        "json_output": "datasets/MySample.json",
        "files": [
            {
                "das_names": [
                    "/DAS/path/to/dataset/NANOAODSIM"
                ],
                "metadata": {
                    "year": "2022_preEE",
                    "isMC": true,
                    "xsec": 1234.0
                }
            }
        ]
    }
}
```

After adding the entry, build the dataset JSON files (inside the singularity with the venv active):

```bash
# Restricting to European sites — recommended from lxplus
pocket-coffea build-datasets --cfg datasets/datasets_definitions.json -o -rs 'T[123]_(FR|IT|DE|BE|CH|UK)_\w+'
```

This generates `datasets/MySample_redirector.json` (and variants) with the actual file lists. Then update `MET_studies_config_lxplus.py` to reference the new dataset JSON and sample name:

```python
datasets={
    "jsons": [f"{localdir}/datasets/MySample_redirector.json"],
    "filter": {"samples": ["MySample"], "year": ["2022_preEE"]},
},
```

For more details on dataset building options (whitelist/blacklist sites, local prefixes, etc.) see the [PocketCoffea datasets documentation](https://pocketcoffea.readthedocs.io/en/latest/datasets.html).
