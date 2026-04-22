#!/bin/bash
#SBATCH --job-name=regroup_by_ADNC
#SBATCH --time=00:10:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1
#SBATCH --output=/data/users/vmuller/0_master_thesis/log/regroup_by_ADNC_%J.out
#SBATCH --error=/data/users/vmuller/0_master_thesis/log/regroup_by_ADNC_%J.err
#SBATCH --partition=pibu_el8

set -euo pipefail

# Define the variables for the directories and required files
source "$SLURM_SUBMIT_DIR/00_config.sh"

# Activate the conda environment from the config file
activate_conda

# Define the output directory for h5ad files regrouped by ADNC status
OUTPUTDIR="$PATH_DONOR_SUBSETS/group"
mkdir -p "$OUTPUTDIR"

# Run the script to prepare the bed files for each cell type
python "$PATH_SCRIPTS/02b_regroup_by_ADNC.py" \
  "$PATH_DONOR_SUBSETS" \
   "$OUTPUTDIR"

echo "Finished"