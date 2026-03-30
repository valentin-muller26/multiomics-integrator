#!/bin/bash
#SBATCH --job-name=download_SEAD
#SBATCH --time=04:00:00
#SBATCH --mem=6G
#SBATCH --cpus-per-task=2
#SBATCH --output=/data/users/vmuller/0_master_thesis/log/download_SEAD_%J.out
#SBATCH --error=/data/users/vmuller/0_master_thesis/log/download_SEAD_%J.err
#SBATCH --partition=pibu_el8

set -euo pipefail

#Define the variables for the directories and required files
source "$SLURM_SUBMIT_DIR/00_config.sh"
mkdir -p $PATH_DATA
cd $PATH_DATA

# RNA
if [ ! -f "SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad" ]; then
    echo "Downloading RNA"
    wget -c https://sea-ad-single-cell-profiling.s3.amazonaws.com/MTG/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad
else
    echo "RNA file already exists, skipping."
fi

# ATAC
if [ ! -f "SEAAD_MTG_ATACseq_final-nuclei.2024-02-13.h5ad" ]; then
    echo "Downloading ATAC"
    wget -c https://sea-ad-single-cell-profiling.s3.amazonaws.com/MTG/ATACseq/previous_objects/SEAAD_MTG_ATACseq_final-nuclei.2024-02-13.h5ad
else
    echo "ATAC file already exists, skipping."
fi

# Download reference 10x GRCh38 genome
mkdir -p $REFDIR
cd $REFDIR
if [ -d "refdata-cellranger-arc-GRCh38-2020-A-2.0.0" ]; then
    echo "The reference genome already exists, skipping."
else
    wget -c https://cf.10xgenomics.com/supp/cell-arc/refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz
    tar -xzf refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz
    #unzip the annotation file
    gunzip -k $REFDIR/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz
    rm refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz
fi