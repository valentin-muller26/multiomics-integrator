#!/bin/bash

export CONDA_ENV=scdesign

#Define the variables for the main directories 
export WORKDIR=/data/users/vmuller/0_master_thesis
export LOGDIR=$WORKDIR/log
export PATH_DATA=$WORKDIR/data/SEAD_Dataset/
export PATH_DONOR_SUBSETS=$WORKDIR/data/SEAD_Dataset/patient_subsets
export PATH_SCRIPTS=$WORKDIR/src/scDesign_pipeline

# Define the URLs and filenames for the RNA and ATAC data from the SEAD dataset.
export RNA_FILE_SEAD="SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad"
export RNA_URL="https://sea-ad-single-cell-profiling.s3.amazonaws.com/MTG/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad"
export ATAC_FILE_SEAD="SEAAD_MTG_ATACseq_final-nuclei.2024-02-13.h5ad"
export ATAC_URL="https://sea-ad-single-cell-profiling.s3.amazonaws.com/MTG/ATACseq/previous_objects/SEAAD_MTG_ATACseq_final-nuclei.2024-02-13.h5ad"


# Precomputed stats for the mapmycells annotation
export PRECOMPUTED_STATS_URL="https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/mapmycells/SEAAD/20240831/precomputed_stats.20231120.sea_ad.MTG.h5"
export PRECOMPUTED_STATS_FILE="$WORKDIR/data/mapmycells/precomputed_stats/precomputed_stats.20231120.sea_ad.MTG.h5"

#Create the directory for the error and output file if not present
mkdir -p $LOGDIR
mkdir -p $PATH_DATA

# Function to activate the conda environment
activate_conda() {
    unset PYTHONPATH
    export CONDA_PREFIX="$HOME/.conda/envs/$CONDA_ENV"
    export PATH="$CONDA_PREFIX/bin:$PATH"
}

