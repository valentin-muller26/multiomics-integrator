#!/bin/bash
#SBATCH --partition=pibu_el8
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --job-name=create_pseudobulk
#SBATCH --output=/data/users/vmuller/0_master_thesis/log/create_pseudobulk_%J.out
#SBATCH --error=/data/users/vmuller/0_master_thesis/log/create_pseudobulk_%J.err

set -euo pipefail

echo $(date)

# Load the configuration variables and activate the conda environment
source "$SLURM_SUBMIT_DIR/00_config.sh"
activate_conda

#Define the input directories for the simulated single-cell RNA and ATAC data and the output directories for the pseudobulk data
INPUT_RNA_DIR="/data/users/vmuller/0_master_thesis/data/SEAD_Dataset/patient_subsets/simulated_data/RNA"
OUTPUT_RNA_DIR="/data/users/vmuller/0_master_thesis/data/SEAD_Dataset/patient_subsets/simulated_pseudobulk/RNA"
INPUT_ATAC_DIR="/data/users/vmuller/0_master_thesis/data/SEAD_Dataset/patient_subsets/simulated_data/ATAC"
OUTPUT_ATAC_DIR="/data/users/vmuller/0_master_thesis/data/SEAD_Dataset/patient_subsets/simulated_pseudobulk/ATAC"


mkdir -p "$OUTPUT_RNA_DIR" "$OUTPUT_ATAC_DIR"

# Find all h5ad files in the input directories
mapfile -t LIST_FILES_RNA < <(find "$INPUT_RNA_DIR" -name "*.h5ad" -print0 | tr '\0' '\n')
mapfile -t LIST_FILES_ATAC < <(find "$INPUT_ATAC_DIR" -name "*.h5ad" -print0 | tr '\0' '\n')

echo "Found ${#LIST_FILES_RNA[@]} RNA files to process."
echo "Found ${#LIST_FILES_ATAC[@]} ATAC files to process."

# Process RNA and ATAC in parallel
for file in "${LIST_FILES_RNA[@]}"; do
    echo "Processing RNA: $file"
    python "$PATH_SCRIPTS/09a_create_pseudobulk.py" \
        --input_path "$file" \
        --output_path "$OUTPUT_RNA_DIR"
done &

for file in "${LIST_FILES_ATAC[@]}"; do
    echo "Processing ATAC: $file"
    python "$PATH_SCRIPTS/09a_create_pseudobulk.py" \
        --input_path "$file" \
        --output_path "$OUTPUT_ATAC_DIR"
done &

wait
echo "Finished"
echo $(date)