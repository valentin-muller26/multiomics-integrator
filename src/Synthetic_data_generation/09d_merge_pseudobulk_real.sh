#!/bin/bash
#SBATCH --partition=pshort_el8
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --job-name=merge_pseudobulk
#SBATCH --output=/data/users/vmuller/0_master_thesis/log/merge_pseudobulk_%J.out
#SBATCH --error=/data/users/vmuller/0_master_thesis/log/merge_pseudobulk_%J.err


set -euo pipefail

echo $(date)

# Load the configuration variables and activate the conda environment
source "$SLURM_SUBMIT_DIR/00_config.sh"
activate_conda

# Define the input directories for the pseudobulk RNA and ATAC data and the output directory for the merged pseudobulk data
INPUT_RNA_DIR="/data/users/vmuller/0_master_thesis/data/SEAD_Dataset/patient_subsets/real_pseudobulk/RNA"
INPUT_ATAC_DIR="/data/users/vmuller/0_master_thesis/data/SEAD_Dataset/patient_subsets/real_pseudobulk/ATAC"
OUTPUT_DIR="/data/users/vmuller/0_master_thesis/data/SEAD_Dataset/patient_subsets/real_pseudobulk_merged"

mkdir -p "$OUTPUT_DIR"

python "$PATH_SCRIPTS/09d_merge_pseudobulk_real.py" \
    "$INPUT_RNA_DIR" \
    "$INPUT_ATAC_DIR" \
    "$OUTPUT_DIR"

echo "Finished"
