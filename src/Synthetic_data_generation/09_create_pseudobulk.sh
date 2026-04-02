#!/bin/bash
#SBATCH --partition=pshort_el8
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --job-name=create_pseudobulk
#SBATCH --output=/data/users/vmuller/0_master_thesis/log/create_pseudobulk_%J.out
#SBATCH --error=/data/users/vmuller/0_master_thesis/log/create_pseudobulk_%J.err

set -euo pipefail

echo $(date)
source "$SLURM_SUBMIT_DIR/00_config.sh"
activate_conda

INPUT_DIR="/data/users/vmuller/0_master_thesis/data/SEAD_Dataset/patient_subsets/simulated_test/RNA"
OUTPUT_DIR="/data/users/vmuller/0_master_thesis/data/SEAD_Dataset/patient_subsets/simulated_test_pseudobulk/RNA"
mkdir -p "$OUTPUT_DIR"

LIST_FILES=($(find "$INPUT_DIR" -name "*.h5ad"))
echo "Found ${#LIST_FILES[@]} files to process."

for file in "${LIST_FILES[@]}"; do
    echo "Processing $file"
    python "$WORKDIR/src/Synthetic_data_generation/create_pseudobulk.py" \
        "$file" \
        "$OUTPUT_DIR"
done



