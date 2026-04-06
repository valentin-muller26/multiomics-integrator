#!/bin/bash
#SBATCH --job-name=cell_type_distribution_analysis
#SBATCH --time=00:10:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=/data/users/vmuller/0_master_thesis/log/cell_type_distribution_analysis_%A_%a.out
#SBATCH --error=/data/users/vmuller/0_master_thesis/log/cell_type_distribution_analysis_%A_%a.err
#SBATCH --partition=pibu_el8

set -euo pipefail

# Define the variables for the directories and required files
source "$SLURM_SUBMIT_DIR/00_config.sh"

# Activate the conda environment from the config file
activate_conda

#Define the output directory for the cell type distribution results
OUTPUTDIR="$PATH_DONOR_SUBSETS/cell_type_distribution"
mkdir -p "$OUTPUTDIR"

# Run the script to analyze the cell type distribution for each patient subset
python "$PATH_SCRIPTS/02a_donor_cell_distribution.py" \
  "$PATH_DONOR_SUBSETS" \
   "$OUTPUTDIR"

echo "Finished"