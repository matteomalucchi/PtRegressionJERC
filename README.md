# PtRegressionJERC

Repository with [PocketCoffea](https://github.com/PocketCoffea/PocketCoffea/tree/main) configurations to compute the MC Truth corrections, JER and MET-Type1 correction for regressed pT jets.

## Setup

### lxplus installation

To setup a local installation on `lxplus`,:

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

# Set the PYTHONPATH to make sure the editable PocketCoffea installation is picked up
export PYTHONPATH=`pwd`

cd ../ 
git clone git@github.com:matteomalucchi/PtRegressionJERC.git
cd PtRegressionJERC
pip install -e ./

# Install the HEPPlotter class
pip install --upgrade  --no-cache-dir git+https://github.com/matteomalucchi/AnalysisConfigs.git
```

After that you should set an alias to activate the PocketCoffea environment because this is called automatically by the `exec.py` script.

On `lxplus`, it can be done by adding the following line to your `.bashrc`:

```bash
alias pocket_coffea='apptainer shell --bind /afs -B /cvmfs/cms.cern.ch \
         --bind /tmp  --bind /eos/cms/ -B /etc/sysconfig/ngbauth-submit \
         -B ${XDG_RUNTIME_DIR}  --env KRB5CCNAME="FILE:${XDG_RUNTIME_DIR}/krb5cc"  \
         /cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-analysis/general/pocketcoffea:lxplus-el9-stable'
```

> [!IMPORTANT]
> For further instructions about the installation of PocketCoffea, you can checkout the [Installation section](https://pocketcoffea.readthedocs.io/en/latest/installation.html#using-apptainer-for-local-development) in the documentation.

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

## Troubleshooting (lxplus)

### Condor jobs fail with `ModuleNotFoundError: No module named 'utils_configs'`

The `configs` package (from [AnalysisConfigs](https://github.com/matteomalucchi/AnalysisConfigs)) is installed by pip into `~/.local/lib/python3.11/site-packages/` instead of the venv, even when the venv is active. This is a known pip behaviour with `--system-site-packages` venvs. Interactive sessions work because Python picks up user site-packages, but condor workers run with `PYTHONNOUSERSITE=1` and cannot find the package.

To verify where `utils_configs` is installed:

```bash
python -c "import utils_configs; print(utils_configs.__file__)"
# Should print: .../pocket_coffea_<env>/lib/python3.11/site-packages/utils_configs/__init__.py
# If it prints ~/.local/..., follow the fix below
```

Fix: remove the user-local install and reinstall into the venv (inside the singularity with the venv active):

```bash
rm -rf ~/.local/lib/python3.11/site-packages/utils_configs
rm -rf ~/.local/lib/python3.11/site-packages/configs-*.dist-info
pip install --no-cache-dir git+https://github.com/matteomalucchi/AnalysisConfigs.git
```

### Changes to PocketCoffea source are not picked up

The singularity image ships its own `pocket_coffea` at `/usr/local/lib/python3.11/site-packages/`, which can take priority over the editable install in the venv. Fix: add a `.pth` file to the venv's site-packages that puts the AFS source directory on `sys.path` before the system install (run once, inside the singularity with the venv active, from the `PocketCoffea` directory):

```bash
echo "$(pwd)" >> $(python3 -c "import site; print(site.getsitepackages()[0])")/pocket_coffea_source.pth
```

Verify the fix:

```bash
python -c "import pocket_coffea; print(pocket_coffea.__file__)"
# Should print: /afs/.../PocketCoffea/pocket_coffea/__init__.py
```
