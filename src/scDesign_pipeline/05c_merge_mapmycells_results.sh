#!/bin/bash
#SBATCH --job-name=merge_mapmycells
#SBATCH --time=04:00:00
#SBATCH --mem=6G
#SBATCH --cpus-per-task=4
#SBATCH --output=/data/users/vmuller/0_master_thesis/log/merge_mapmycells_%J.out
#SBATCH --error=/data/users/vmuller/0_master_thesis/log/merge_mapmycells_%J.err
#SBATCH --partition=pibu_el8

set -euo pipefail

echo $(date)

#Define the variables for the directories and required files
source "$SLURM_SUBMIT_DIR/00_config.sh"

#Activate the conda environment from the config file
activate_conda

#Define the output directory for the merged mapmycells results
MAIN_INPUT_DIR="$PATH_DONOR_SUBSETS/mapmycells_results"
OUTPUTDIR="$PATH_DONOR_SUBSETS/mapmycells_results/merged_results"
mkdir -p "$OUTPUTDIR"

ADNC_STATUS=("Not AD" "Low" "Intermediate" "High")
# Loop through the ADNC groups and merge the results
for STATUS in "${ADNC_STATUS[@]}"; do
    INPUT_DIR="$MAIN_INPUT_DIR/$STATUS/graphs"
    if [[ -d "$INPUT_DIR" ]]; then
        echo "Merging results for ADNC group: $STATUS"
        python "$PATH_SCRIPTS/05c_merge_mapmycells_results.py" \
            --input "$INPUT_DIR" \
            --output "$OUTPUTDIR"
    else
        echo "Warning: Input directory not found for ADNC group $STATUS: $INPUT_DIR"
    fi
done    
