# MET Studies for PNet regressed pT jets

Repository to compute MET Type-1 corrections for PNet regressed pT jets, structured as an analysis configurations for a specific development branch from a fork [PocketCoffea](https://github.com/PocketCoffea/PocketCoffea/tree/main).

## Workflow

### Running the analysis

To run the analysis on Tier3, use the following command:

```bash
pocket-coffea run --cfg MET_studies_config.py --custom-run-options params/t3_run_options_big.yaml -o <output-dir> -e dask@T3_CH_PSI
```

To produce the response plots, use:

```bash
submit_job_10min_25gb_8cpu python plot_MET.py -i <input-dir> -w 8 --histo --novars -o <output-plot-dir>
```
