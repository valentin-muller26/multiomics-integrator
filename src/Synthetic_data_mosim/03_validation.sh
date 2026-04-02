#!/bin/bash
#SBATCH --job-name=mosim_validation
#SBATCH --time=1:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --output=/data/users/vmuller/0_master_thesis/log/mosim_validation_%J.out
#SBATCH --error=/data/users/vmuller/0_master_thesis/log/mosim_validation_%J.err
#SBATCH --partition=pshort_el8

set -euo pipefail
echo $(date)
source "$SLURM_SUBMIT_DIR/00_config.sh"

activate_conda

DONOR_ID="HC19"
INPUTDIR="/data/users/vmuller/0_master_thesis/data/data_interleukines/simulated_data"
OUTPUTDIR="/data/users/vmuller/0_master_thesis/data/data_interleukines/validation_results"

mkdir -p "$OUTPUTDIR"

Rscript "$WORKDIR/src/Synthetic_data_mosim/validation_mosim.R" \
    --donor_name "$DONOR_ID" \
    --inputdir "$INPUTDIR" \
    --outdir "$OUTPUTDIR"