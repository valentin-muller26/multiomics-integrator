#!/bin/bash
export WORKDIR=/data/users/vmuller/0_master_thesis
export LOGDIR=$WORKDIR/log
export CONDA_ENV=multiomic_analysis
export PATH_DATA=$WORKDIR/data/data_interleukines
export PATH_SCRIPTS=$WORKDIR/src/MultiOmic_analysis

#Create the directory for the error and output file if not present
mkdir -p $LOGDIR

# Function to activate the conda environment
activate_conda() {
    unset PYTHONPATH
    export CONDA_PREFIX="/home/vmuller/.conda/envs/$CONDA_ENV"
    export PATH="$CONDA_PREFIX/bin:$PATH"
}