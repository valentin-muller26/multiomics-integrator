#!/bin/bash
#SBATCH --job-name=extracting_patient
#SBATCH --time=12:00:00
#SBATCH --mem=120G
#SBATCH --cpus-per-task=1
#SBATCH --output=/data/users/vmuller/0_master_thesis/log/extracting_patient_%J.out
#SBATCH --error=/data/users/vmuller/0_master_thesis/log/extracting_patient_%J.err
#SBATCH --partition=pibu_el8

set -euo pipefail

#Define the variables for the directories and required files
source "$SLURM_SUBMIT_DIR/00_config.sh"

#Activate the conda environment from the config file
activate_conda

RNAFILENAME="SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad"
ATACFILENAME="SEAAD_MTG_ATACseq_final-nuclei.2024-02-13.h5ad"
OUTPUTDIR="$PATH_DATA/patient_subsets"

# Run the script
python $WORKDIR/src/Synthetic_data_generation/extracting_patient.py  $PATH_DATA/$RNAFILENAME $PATH_DATA/$ATACFILENAME $OUTPUTDIR

echo "Finished"