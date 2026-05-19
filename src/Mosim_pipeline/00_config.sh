#!/bin/bash
export WORKDIR=/data/users/vmuller/0_master_thesis
export LOGDIR=$WORKDIR/log
export CONDA_ENV=mosim
export PATH_DATA=$WORKDIR/data/data_interleukines
export PATH_SCRIPTS=$WORKDIR/src/Mosim_pipeline
#Create the directory for the error and output file if not present
mkdir -p "$LOGDIR"
mkdir -p "$PATH_DATA"

# Function to activate the conda environment
activate_conda() {
    unset PYTHONPATH
    export CONDA_PREFIX="$HOME/.conda/envs/$CONDA_ENV"
    export PATH="$CONDA_PREFIX/bin:$PATH"
}

