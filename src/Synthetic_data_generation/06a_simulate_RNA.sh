#!/bin/bash
#SBATCH --job-name=sim_scDesignRNA
#SBATCH --time=4-00:00:00
#SBATCH --mem=80G
#SBATCH --cpus-per-task=3
#SBATCH --output=/data/users/vmuller/0_master_thesis/log/sim_scDesignRNA_%A_%a.out
#SBATCH --error=/data/users/vmuller/0_master_thesis/log/sim_scDesignRNA_%A_%a.err
#SBATCH --partition=pibu_el8
#SBATCH --array=0-3
set -euo pipefail
 
echo $(date)
source "$SLURM_SUBMIT_DIR/00_config.sh"
activate_conda
export RETICULATE_PYTHON="$CONDA_PREFIX/bin/python"

MANIFEST="$WORKDIR/data/SEAD_Dataset/models/RNA/manifest.csv"
MODEL_DIR=$(awk -v line="$((SLURM_ARRAY_TASK_ID + 1))" -F',' 'NR==line {print $2}' "$MANIFEST")
GROUPNAME=$(awk -v line="$((SLURM_ARRAY_TASK_ID + 1))" -F',' 'NR==line {print $1}' "$MANIFEST")
SIM_OUT="$WORKDIR/data/SEAD_Dataset/patient_subsets/simulated_data/RNA/$GROUPNAME"
 
mkdir -p "$SIM_OUT"
 
Rscript "$WORKDIR/src/Synthetic_data_generation/simulate_RNA.R" \
  --modeldir  "$MODEL_DIR" \
  --outdir    "$SIM_OUT" \
  --donor_id  "$GROUPNAME"\
  --n_rep 50
 
echo "Finished"
echo $(date)

