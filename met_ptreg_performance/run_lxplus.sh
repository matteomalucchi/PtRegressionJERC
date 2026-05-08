#!/bin/bash
set -e

OUTPUT_DIR=${1:?Usage: $0 <output-dir>}

USER_INITIAL=${USER:0:1}
AFS_USER=/afs/cern.ch/user/${USER_INITIAL}/${USER}
EOS_USER=/eos/user/${USER_INITIAL}/${USER}

WORK_DIR=${AFS_USER}/PtRegressionJERC/met_ptreg_performance
VENV=${AFS_USER}/PocketCoffea/pocket_coffea_env/bin/activate
PYTHONPATH_DIR=${AFS_USER}/PocketCoffea
# IMAGE=/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-analysis/general/pocketcoffea:lxplus-el9-stable

source ${VENV}
export PYTHONPATH=${PYTHONPATH_DIR}:\$PYTHONPATH
cd ${WORK_DIR}
pocket-coffea run \
    --cfg MET_studies_config_lxplus.py \
    --custom-run-options ./params/lxplus_run_options_small.yaml \
    -o ${OUTPUT_DIR} \
    -e condor@lxplus


# apptainer exec \
#     --bind /afs \
#     -B /cvmfs/ \
#     -B ${EOS_USER} \
#     --bind /tmp \
#     --bind /eos/cms/ \
#     -B /etc/sysconfig/ngbauth-submit \
#     -B ${XDG_RUNTIME_DIR} \
#     --env KRB5CCNAME="FILE:${XDG_RUNTIME_DIR}/krb5cc" \
#     "$IMAGE" \
#     bash -c "
#     "
