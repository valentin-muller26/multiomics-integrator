#!/bin/bash
#SBATCH --job-name=download_SEAD_dataset
#SBATCH --time=04:00:00
#SBATCH --mem=6G
#SBATCH --cpus-per-task=2
#SBATCH --output=/data/users/vmuller/0_master_thesis/log/download_SEAD_dataset_%J.out
#SBATCH --error=/data/users/vmuller/0_master_thesis/log/download_SEAD_dataset_%J.err
#SBATCH --partition=pibu_el8

set -euo pipefail


source "$SLURM_SUBMIT_DIR/00_config.sh"

#Define the variables for the output directories
mkdir -p $PATH_DATA
cd $PATH_DATA


# Download RNA donor data 
if [ ! -f "$RNA_FILE_SEAD" ]; then
    echo "Downloading RNA"
    wget -c "$RNA_URL"
else
    echo "RNA file already exists, skipping."
fi

# Download ATAC donor data
if [ ! -f "$ATAC_FILE_SEAD" ]; then
    echo "Downloading ATAC"
    wget -c "$ATAC_URL"
else
    echo "ATAC file already exists, skipping."
fi

echo "Finished"