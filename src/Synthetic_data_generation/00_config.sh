#!/bin/bash

export CONDA_ENV=scReadSim

#Define the variables for the main directories 
export WORKDIR=/data/users/vmuller/0_master_thesis
export LOGDIR=$WORKDIR/log
export PATH_DATA=$WORKDIR/data/SEAD_Dataset/
export PATH_DONOR_SUBSETS=$WORKDIR/data/SEAD_Dataset/patient_subsets

# Define the URLs and filenames for the RNA and ATAC data from the SEAD dataset.
export RNA_FILE_SEAD="SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad"
export RNA_URL="https://sea-ad-single-cell-profiling.s3.amazonaws.com/MTG/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad"
export ATAC_FILE_SEAD="SEAAD_MTG_ATACseq_final-nuclei.2024-02-13.h5ad"
export ATAC_URL="https://sea-ad-single-cell-profiling.s3.amazonaws.com/MTG/ATACseq/previous_objects/SEAAD_MTG_ATACseq_final-nuclei.2024-02-13.h5ad"

#Create the directory for the error and output file if not present
mkdir -p $LOGDIR

# Function to activate the conda environment
activate_conda() {
    unset PYTHONPATH
    export CONDA_PREFIX="/data/users/vmuller/miniforge3/envs/$CONDA_ENV"
    export PATH="$CONDA_PREFIX/bin:$PATH"
}

