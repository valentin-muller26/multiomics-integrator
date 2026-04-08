#!/bin/bash
#SBATCH --job-name=fit_scDesign3_ATAC
#SBATCH --time=3-00:00:00
#SBATCH --mem=120G
#SBATCH --cpus-per-task=2
#SBATCH --output=/data/users/vmuller/0_master_thesis/log/fit_scDesignATAC_%A_%a.out
#SBATCH --error=/data/users/vmuller/0_master_thesis/log/fit_scDesignATAC_%A_%a.err
#SBATCH --partition=pibu_el8
#SBATCH --array=0-3

set -euo pipefail
 
echo $(date)

# Load the configuration and activate the conda environment
source "$SLURM_SUBMIT_DIR/00_config.sh"
activate_conda
 
# Retrieve the donor name based on the SLURM_ARRAY_TASK_ID
PATH_DONOR_DATA="$WORKDIR/data/SEAD_Dataset/patient_subsets"

# Get the ATAC file and group name for the current task ID
GROUP="$PATH_DONOR_DATA/group/ATAC/manifest.csv"
ATAC_FILE=$(awk -v line="$((SLURM_ARRAY_TASK_ID + 1))" -F',' 'NR==line {print $2}' "$GROUP")
GROUPNAME=$(awk -v line="$((SLURM_ARRAY_TASK_ID + 1))" -F',' 'NR==line {print $1}' "$GROUP")
MODEL_DIR="$WORKDIR/data/SEAD_Dataset/models/ATAC/$GROUPNAME"

echo "$ATAC_FILE"
echo "$MODEL_DIR"
mkdir -p "$MODEL_DIR"
 
#Run the R script to fit the scDesign3 model for the current donor subset
Rscript "$PATH_SCRIPTS/fit_model_scDesignAtac.R" \
  --atacfile "$ATAC_FILE" \
  --outdir  "$MODEL_DIR"


# Update the manifest file with the new model information
# Use a lock file to prevent concurrent writes to the manifest
MANIFEST="$WORKDIR/data/SEAD_Dataset/models/ATAC/manifest.csv"
LOCKFILE="$MANIFEST.lock"

(
    flock -x 200
    if ! grep -q "^$GROUPNAME," "$MANIFEST" 2>/dev/null; then
        echo "$GROUPNAME,$MODEL_DIR" >> "$MANIFEST"
    fi
) 200>"$LOCKFILE"