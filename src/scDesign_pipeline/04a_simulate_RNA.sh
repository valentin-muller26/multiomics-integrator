#!/bin/bash
#SBATCH --job-name=sim_scDesignRNA
#SBATCH --time=8-00:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=5
#SBATCH --output=/data/users/vmuller/0_master_thesis/log/sim_scDesignRNA_%A_%a.out
#SBATCH --error=/data/users/vmuller/0_master_thesis/log/sim_scDesignRNA_%A_%a.err
#SBATCH --partition=pibu_el8
#SBATCH --array=0-3


set -euo pipefail 
echo "$(date)"

# Load config and activate conda environment
source "$SLURM_SUBMIT_DIR/00_config.sh"
activate_conda
export RETICULATE_PYTHON="$CONDA_PREFIX/bin/python"

# Load the manifest indicating which model to use for each group
MANIFEST="$PATH_DATA/models/RNA/manifest.csv"

# Extract the model directory and group name for the current task
MODEL_DIR=$(awk -v line="$((SLURM_ARRAY_TASK_ID + 1))" -F',' 'NR==line {print $2}' "$MANIFEST")
GROUPNAME=$(awk -v line="$((SLURM_ARRAY_TASK_ID + 1))" -F',' 'NR==line {print $1}' "$MANIFEST")

# Define the output directory for the simulated data
SIM_OUT="$PATH_DONOR_SUBSETS/simulated_data/RNA/$GROUPNAME"
mkdir -p "$SIM_OUT"
 
# Run the R script to simulate RNA data using scDesign3
Rscript "$PATH_SCRIPTS/04a_simulate_RNA.R" \
  --modeldir  "$MODEL_DIR" \
  --outdir    "$SIM_OUT" \
  --ADNC_cat  "$GROUPNAME"\
  --n_rep 1
  
 
echo "Finished"
echo $(date)

