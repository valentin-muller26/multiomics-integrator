#!/bin/bash
#SBATCH --job-name=fit_scDesignRNA
#SBATCH --time=1-00:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=2
#SBATCH --output=/data/users/vmuller/0_master_thesis/log/fit_scDesignRNA_%A_%a.out
#SBATCH --error=/data/users/vmuller/0_master_thesis/log/fit_scDesignRNA_%A_%a.err
#SBATCH --partition=pibu_el8
#SBATCH --array=0-3
set -euo pipefail
 
echo $(date)
source "$SLURM_SUBMIT_DIR/00_config.sh"
activate_conda
 

# Retrieve the donor name based on the SLURM_ARRAY_TASK_ID
PATH_DONOR_DATA="$WORKDIR/data/SEAD_Dataset/patient_subsets"

GROUP="$PATH_DONOR_DATA/group/RNA/manifest.csv"
RNA_FILE=$(awk -v line="$((SLURM_ARRAY_TASK_ID + 1))" -F',' 'NR==line {print $2}' "$GROUP")
GROUPNAME=$(awk -v line="$((SLURM_ARRAY_TASK_ID + 1))" -F',' 'NR==line {print $1}' "$GROUP")
MODEL_DIR="$WORKDIR/data/SEAD_Dataset/models/RNA/$GROUPNAME"

echo "$RNA_FILE"
echo "$MODEL_DIR"
mkdir -p "$MODEL_DIR"
 
Rscript "$WORKDIR/src/Synthetic_data_generation/fit_model_scDesignRNA.R" \
  --rnafile "$RNA_FILE" \
  --outdir  "$MODEL_DIR"

MANIFEST="$WORKDIR/data/SEAD_Dataset/models/RNA/manifest.csv"
LOCKFILE="$MANIFEST.lock"

(
    flock -x 200
    if ! grep -q "^$GROUPNAME," "$MANIFEST" 2>/dev/null; then
        echo "$GROUPNAME,$MODEL_DIR" >> "$MANIFEST"
    fi
) 200>"$LOCKFILE"