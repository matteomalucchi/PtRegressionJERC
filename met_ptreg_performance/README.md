# MET Studies for regressed pT jets

Repository to compute the response and resolution of MET as well as the effect of the Type-1 correction derived for PNet regressed pT jets.

This repository is structured as an analysis configurations for [PocketCoffea](https://github.com/PocketCoffea/PocketCoffea/tree/main) and was originally hosted [here](https://github.com/matteomalucchi/AnalysisConfigs/tree/legacy-met-config/configs/MET_studies).

## Workflow

In order to reduce the 


### Running the analysis

To run the analysis on Tier3, use the following command:

```bash
pocket-coffea run --cfg MET_studies_config.py --custom-run-options params/t3_run_options_big.yaml -o <output-dir> -e dask@T3_CH_PSI
```

To produce the response plots, use:

```bash
submit_job_10min_25gb_8cpu python plot_MET.py -i <input-dir> -w 8 --histo --novars -o <output-plot-dir>
```

### Lxplus

To run on lxplus, first activate the singularity and the venv. After you can pocket coffea using the `lxplus` executors `condor@lxplus` or `dask@lxplus`:

```bash
# activate singularity
pocket_coffea

# activate venv
cd PocketCoffea
source pocket_coffea_env/bin/activate
export PYTHONPATH=`pwd`

# run pocket coffea
cd ../PtRegressionJERC/met_ptreg_performance 
pocket-coffea run --cfg MET_studies_config_lxplus.py --custom-run-options ./params/lxplus_run_options_big.yaml -o <output-dir> -e condor@lxplus
```
