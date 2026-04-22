#!/bin/bash
#SBATCH --partition=pibu_el8
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=300G
#SBATCH --job-name=validate
#SBATCH --output=/data/users/vmuller/0_master_thesis/log/validate_simulation_%J.out
#SBATCH --error=/data/users/vmuller/0_master_thesis/log/validate_simulation_%J.err

set -euo pipefail
 
echo $(date)
source "$SLURM_SUBMIT_DIR/00_config.sh"
activate_conda




#RNA_FILE="/data/users/vmuller/0_master_thesis/data/SEAD_Dataset/patient_subsets/group/RNA/Not AD.h5ad"
RNA_FILE="/data/users/vmuller/0_master_thesis/data/SEAD_Dataset/patient_subsets/scDesign3_models/NOTAD_RNA/training_cells.h5ad"
SIMULATE_FILE="/data/users/vmuller/0_master_thesis/data/SEAD_Dataset/patient_subsets/simulated_test/NOT_AD_rep01_seed101_scDesign3.h5ad"
OUTPUTDIR="/data/users/vmuller/0_master_thesis/data/SEAD_Dataset/patient_subsets/validation/RNA"
 
mkdir -p "$OUTPUTDIR"

Rscript  "$PATH_SCRIPTS/08_validate_simulation.R" \
  --real_sce      "$RNA_FILE" \
  --sim_sce       "$SIMULATE_FILE" \
  --modality      RNA \
  --output_dir    "$OUTPUTDIR" \
  --prefix        ADNC_NOTAD \
  --cell_type_col cell_type \
  --n_cores       8
