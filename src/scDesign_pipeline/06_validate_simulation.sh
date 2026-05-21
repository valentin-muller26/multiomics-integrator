#!/bin/bash
#SBATCH --partition=pibu_el8
#SBATCH --time=4-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=300G
#SBATCH --job-name=validate
#SBATCH --output=/data/users/vmuller/0_master_thesis/log/validate_%A_%a.out
#SBATCH --error=/data/users/vmuller/0_master_thesis/log/validate_%A_%a.err
#SBATCH --array=0-7

echo "Starting validation for array index: $SLURM_ARRAY_TASK_ID"
echo "$(date)"

# Load configuration and activate conda environment
set -euo pipefail
source "$SLURM_SUBMIT_DIR/00_config.sh"
activate_conda

# --- Configuration ---
MODALITIES=("RNA" "ATAC")
ADNC_STATUS=("Not AD" "Low" "Intermediate" "High")
N_REPS=5

SIM_BASE="$PATH_DONOR_SUBSETS/simulated_data"
OUT_BASE="$PATH_DONOR_SUBSETS/validation"

# Retrieve modality and ADNC status based on the array index
modality_index=$(( SLURM_ARRAY_TASK_ID / 4 ))
adnc_index=$(( SLURM_ARRAY_TASK_ID % 4 ))
MODALITY="${MODALITIES[$modality_index]}"
ADNC="${ADNC_STATUS[$adnc_index]}"

# Define the directory containing the simulated files for the current modality and ADNC status
SIM_DIR="$SIM_BASE/${MODALITY}/${ADNC}"

# List and sort the simulated files, then select the top N_REPS files
mapfile -t SIM_FILES < <(
    find "$SIM_DIR" -maxdepth 1 -name "*_scDesign3.h5ad" -type f \
    | sort -V \
    | head -n "$N_REPS"
)


MODEL_FILE="$PATH_DATA/models/${MODALITY}/$ADNC/training_cells.h5ad"

# Validate each simulated file against the reference model
echo "Validating simulated files for modality=${MODALITY}, ADNC=${ADNC}"
for SIM_FILE in "${SIM_FILES[@]}"; do
    echo "    Simulated file: $SIM_FILE"
    echo "    Reference model: $MODEL_FILE"
    rep_tag=$(echo "$SIM_FILE" | grep -oE 'rep[0-9]+')
    OUTPUTDIR="$OUT_BASE/${MODALITY}/${ADNC}/${rep_tag}"
    mkdir -p "$OUTPUTDIR"

    Rscript "$PATH_SCRIPTS/06_validate_simulation.R" \
        --sim_sce "$SIM_FILE" \
        --real_sce "$MODEL_FILE" \
        --modality "$MODALITY" \
        --prefix "$MODALITY"_"$ADNC"_"$rep_tag" \
        --output_dir "$OUTPUTDIR"
done

echo "Validation completed for modality=${MODALITY}, ADNC=${ADNC}"
echo "$(date)"