#!/bin/bash
#SBATCH --job-name=mosim_validation
#SBATCH --time=1:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --output=/data/users/vmuller/0_master_thesis/log/mosim_validation_%J.out
#SBATCH --error=/data/users/vmuller/0_master_thesis/log/mosim_validation_%J.err
#SBATCH --partition=pshort_el8

# Load the configuration and activate the conda environment
set -euo pipefail
echo "$(date)"
source "$SLURM_SUBMIT_DIR/00_config.sh"
activate_conda

DONOR_ID="HC19"
INPUTDIR="$PATH_DATA/simulated_data"
OUTPUTDIR="$PATH_DATA/validation_results"

mkdir -p "$OUTPUTDIR"

Rscript "$PATH_SCRIPTS/03_validation_mosim.R" \
    --donor_ID "$DONOR_ID" \
    --inputdir "$INPUTDIR" \
    --outdir "$OUTPUTDIR"

echo "Validation completed"
echo "$(date)"