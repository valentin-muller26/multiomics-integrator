#!/bin/bash
#SBATCH --job-name=extract_donor
#SBATCH --time=12:00:00
#SBATCH --mem=120G
#SBATCH --cpus-per-task=1
#SBATCH --output=/data/users/vmuller/0_master_thesis/log/extract_donor_%J.out
#SBATCH --error=/data/users/vmuller/0_master_thesis/log/extract_donor_%J.err
#SBATCH --partition=pibu_el8

set -euo pipefail

#Define the variables for the directories and required files
source "$SLURM_SUBMIT_DIR/00_config.sh"

#Activate the conda environment from the config file
activate_conda


#Define the filenames for the RNA and ATAC data from the SEAD dataset and the output directory for the extracted patient subsets.
OUTPUTDIR="$PATH_DONOR_SUBSETS"

# Run the script
python $PATH_SCRIPTS/02_extract_donor.py  $PATH_DATA/$RNA_FILE_SEAD $PATH_DATA/$ATAC_FILE_SEAD $OUTPUTDIR

echo "Finished"