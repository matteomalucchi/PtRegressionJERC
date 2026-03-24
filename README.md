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
