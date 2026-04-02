#!/bin/bash
export WORKDIR=/data/users/vmuller/0_master_thesis
export APPTAINER_SAMTOOLS_PATH=/containers/apptainer/samtools-1.19.sif

export LOGDIR=$WORKDIR/log
export REFDIR=$WORKDIR/data/refgenome_dir
export CONDA_ENV=scReadSim
export PATH_DATA=$WORKDIR/data/SEAD_Dataset/
export GTF=$WORKDIR/data/refgenome_dir/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf
#Create the directory for the error and output file if not present
mkdir -p $LOGDIR

# Function to activate the conda environment
activate_conda() {
    unset PYTHONPATH
    export CONDA_PREFIX="/data/users/vmuller/miniforge3/envs/$CONDA_ENV"
    export PATH="$CONDA_PREFIX/bin:$PATH"
}

