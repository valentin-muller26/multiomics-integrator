#!/bin/bash
#SBATCH --job-name=mosim_simulation
#SBATCH --time=1:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --output=/data/users/vmuller/0_master_thesis/log/mosim_simulation_%J.out
#SBATCH --error=/data/users/vmuller/0_master_thesis/log/mosim_simulation_%J.err
#SBATCH --partition=pshort_el8

set -euo pipefail
echo "$(date)"
source "$SLURM_SUBMIT_DIR/00_config.sh"
activate_conda


DONOR_ID="DA17"
OUTPUTDIR="$PATH_DATA/simulated_data"

mkdir -p "$OUTPUTDIR"

Rscript "$PATH_SCRIPTS/02_run_mosim.R" \
    --donor_ID "$DONOR_ID" \
    --inputdir  "$PATH_DATA" \
    --outdir "$OUTPUTDIR"

echo "Simulation completed"
echo "$(date)"