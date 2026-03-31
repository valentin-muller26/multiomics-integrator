#!/bin/bash
#SBATCH --job-name=regroup_by_ADNC
#SBATCH --time=00:10:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1
#SBATCH --output=/data/users/vmuller/0_master_thesis/log/regroup_by_ADNC_%A_%a.out
#SBATCH --error=/data/users/vmuller/0_master_thesis/log/regroup_by_ADNC_%A_%a.err
#SBATCH --partition=pibu_el8

set -euo pipefail

# Define the variables for the directories and required files
source "$SLURM_SUBMIT_DIR/00_config.sh"

# Activate the conda environment from the config file
activate_conda

PATH_DONOR_DATA="$WORKDIR/data/SEAD_Dataset/patient_subsets"
OUTPUT_DIR="$WORKDIR/data/SEAD_Dataset/patient_subsets/group"

# Run the script to prepare the bed files for each cell type
python "$WORKDIR/src/Synthetic_data_generation/regroup_by_ADNC.py" \
  "$PATH_DONOR_DATA" \
   "$OUTPUT_DIR"