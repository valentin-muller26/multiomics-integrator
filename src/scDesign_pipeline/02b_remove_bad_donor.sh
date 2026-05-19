#!/bin/bash
#SBATCH --job-name=remove_bad_donor
#SBATCH --time=00:10:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=/data/users/vmuller/0_master_thesis/log/remove_bad_donor_%J.out
#SBATCH --error=/data/users/vmuller/0_master_thesis/log/remove_bad_donor_%J.err
#SBATCH --partition=pibu_el8

set -euo pipefail
echo "$(date)"
# Define the variables for the directories and required files
source "$SLURM_SUBMIT_DIR/00_config.sh"
# Activate the conda environment from the config file
activate_conda

# Define the output directory for h5ad files with bad donors removed
OUTPUTDIR="$PATH_DONOR_SUBSETS/bad_donor_removed"
mkdir -p "$OUTPUTDIR"

FILE_BAD_DONORS="$PATH_DATA/bad_donors.txt"

while IFS= read -r DONOR || [[ -n "$DONOR" ]]; do
  echo "Processing donor: $DONOR"
   RNA_FILE="$PATH_DONOR_SUBSETS/${DONOR}_RNA_multiome_subset.h5ad"
    ATAC_FILE="$PATH_DONOR_SUBSETS/${DONOR}_ATAC_multiome_subset.h5ad"
  if [[ -f "$RNA_FILE" ]]; then
    mv "$RNA_FILE" "$OUTPUTDIR"
    echo "Moved $RNA_FILE to $OUTPUTDIR"
  else
    echo "Warning: RNA file $RNA_FILE does not exist. Skipping donor $DONOR."
  fi
  if [[ -f "$ATAC_FILE" ]]; then
    mv "$ATAC_FILE" "$OUTPUTDIR"
    echo "Moved $ATAC_FILE to $OUTPUTDIR"
  else
    echo "Warning: ATAC file $ATAC_FILE does not exist. Skipping donor $DONOR."
  fi
done < "$FILE_BAD_DONORS"

echo "Finished removing bad donors"
echo "$(date)"