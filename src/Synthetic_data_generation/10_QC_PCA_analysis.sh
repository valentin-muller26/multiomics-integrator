#!/bin/bash
#SBATCH --partition=pibu_el8
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --job-name=pca_real_vs_sim
#SBATCH --output=/data/users/vmuller/0_master_thesis/log/pca_real_vs_sim_%J.out
#SBATCH --error=/data/users/vmuller/0_master_thesis/log/pca_real_vs_sim_%J.err

set -euo pipefail

echo $(date)
# Load config and activate conda environment
source "$SLURM_SUBMIT_DIR/00_config.sh"
activate_conda

#Define the paths to the real and simulated pseudobulk data
REAL_RNA="$PATH_DONOR_SUBSETS/real_pseudobulk_merged/real_RNA_merged_pseudobulk.txt"
REAL_ATAC="$PATH_DONOR_SUBSETS/real_pseudobulk_merged/real_ATAC_merged_pseudobulk.txt"
SIM_RNA="$PATH_DONOR_SUBSETS/simulated_pseudobulk_merged/RNA_merged_pseudobulk.txt"
SIM_ATAC="$PATH_DONOR_SUBSETS/simulated_pseudobulk_merged/ATAC_merged_pseudobulk.txt"

# Define the output directory for the PCA results
OUTPUT="$PATH_DONOR_SUBSETS/PCA"
mkdir -p "$OUTPUT"

# Run the R script to perform PCA analysis comparing real and simulated data
Rscript "$PATH_SCRIPTS/10_QC_PCA_analysis.R" \
    --realRNA "$REAL_RNA" \
    --simRNA "$SIM_RNA" \
    --realATAC "$REAL_ATAC" \
    --simATAC "$SIM_ATAC" \
    --outdir "$OUTPUT"

echo "Finished"
echo $(date)