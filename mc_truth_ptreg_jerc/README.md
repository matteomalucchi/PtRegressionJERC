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

## Contributors

- [Matteo Malucchi](https://github.com/matteomalucchi)
- [Jessy Daniel](https://github.com/jessy-daniel)
- [Jason Guo](https://github.com/erfz)
