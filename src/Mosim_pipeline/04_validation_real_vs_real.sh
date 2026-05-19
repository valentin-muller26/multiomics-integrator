#!/bin/bash
#SBATCH --job-name=mosim_validation_real_vs_real
#SBATCH --time=1:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --output=/data/users/vmuller/0_master_thesis/log/mosim_validation_real_vs_real_%J.out
#SBATCH --error=/data/users/vmuller/0_master_thesis/log/mosim_validation_real_vs_real_%J.err
#SBATCH --partition=pshort_el8

# Load the configuration and activate the conda environment
set -euo pipefail
echo "$(date)"
source "$SLURM_SUBMIT_DIR/00_config.sh"
activate_conda

INPUTDIR="$PATH_DATA/simulated_data/real"
DONOR1="HC19"
DONOR2="HC17"
OUTPUTDIR="$PATH_DATA/validation_results/real_vs_real"
mkdir -p "$OUTPUTDIR"

Rscript "$PATH_SCRIPTS/04_validation_real_vs_real.R" \
    --donor_1 "$DONOR1" \
    --donor_2 "$DONOR2" \
    --input_path "$INPUTDIR" \
    --outdir "$OUTPUTDIR" \
    --modality "ATAC"

echo "Validation ATAC completed"

Rscript "$PATH_SCRIPTS/04_validation_real_vs_real.R" \
    --donor_1 "$DONOR1" \
    --donor_2 "$DONOR2" \
    --input_path "$INPUTDIR" \
    --outdir "$OUTPUTDIR" \
    --modality "RNA"

echo "Validation RNA completed"
echo "$(date)"