#!/bin/bash
#SBATCH --partition=pibu_el8
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --job-name=create_pseudobulkreal
#SBATCH --output=/data/users/vmuller/0_master_thesis/log/create_pseudobulkreal_%J.out
#SBATCH --error=/data/users/vmuller/0_master_thesis/log/create_pseudobulkreal_%J.err

set -euo pipefail

echo $(date)

# Load the configuration variables and activate the conda environment
source "$SLURM_SUBMIT_DIR/00_config.sh"
activate_conda


#Retrieve the training data in the models folder and rename them to contain the group name, 
#and copy to a new folder corresponding to there omic type (RNA or ATAC) in the training data directory
TRAINING_DATA_DIR="$PATH_DATA/training"
MODEL_DIR="$WORKDIR/data/SEAD_Dataset/models"
mkdir -p "$TRAINING_DATA_DIR"
mkdir -p "$TRAINING_DATA_DIR/ATAC"
mkdir -p "$TRAINING_DATA_DIR/RNA"
for omic in ATAC RNA; do
    find "$MODEL_DIR/$omic" -name "training_cells.h5ad" | while read -r file; do
        group=$(basename "$(dirname "$file")")
        cp "$file" "$TRAINING_DATA_DIR/$omic/real_${group}.h5ad"
    done
done


#Define the input directories for the real single-cell RNA and ATAC data and the output directories for the pseudobulk data
INPUT_RNA_DIR="$TRAINING_DATA_DIR/RNA"
OUTPUT_RNA_DIR="$PATH_DONOR_SUBSETS/real_pseudobulk/RNA"
INPUT_ATAC_DIR="$TRAINING_DATA_DIR/ATAC"
OUTPUT_ATAC_DIR="$PATH_DONOR_SUBSETS/real_pseudobulk/ATAC"
mkdir -p "$OUTPUT_RNA_DIR" "$OUTPUT_ATAC_DIR"

# Find all h5ad files in the input directories
mapfile -t LIST_FILES_RNA < <(find "$INPUT_RNA_DIR" -name "*.h5ad" -print0 | tr '\0' '\n')
mapfile -t LIST_FILES_ATAC < <(find "$INPUT_ATAC_DIR" -name "*.h5ad" -print0 | tr '\0' '\n')
echo "Found ${#LIST_FILES_RNA[@]} RNA files to process."
echo "Found ${#LIST_FILES_ATAC[@]} ATAC files to process."

# Process RNA and ATAC in parallel
for file in "${LIST_FILES_RNA[@]}"; do
    echo "Processing RNA: $file"
    python "$PATH_SCRIPTS/09c_create_pseudobulk_real.py" \
        --input_path "$file" \
        --output_path "$OUTPUT_RNA_DIR"\
        --modality RNA
done &

for file in "${LIST_FILES_ATAC[@]}"; do
    echo "Processing ATAC: $file"
    python "$PATH_SCRIPTS/09c_create_pseudobulk_real.py" \
        --input_path "$file" \
        --output_path "$OUTPUT_ATAC_DIR"\
        --modality ATAC
done &

wait
echo "Finished"
echo $(date)