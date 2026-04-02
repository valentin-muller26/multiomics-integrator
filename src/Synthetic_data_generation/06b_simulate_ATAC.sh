#!/bin/bash
#SBATCH --job-name=sim_ATAC
#SBATCH --time=1-00:00:00
#SBATCH --mem=80G
#SBATCH --cpus-per-task=2
#SBATCH --output=/data/users/vmuller/0_master_thesis/log/sim_ATAC_%A_%a.out
#SBATCH --error=/data/users/vmuller/0_master_thesis/log/sim_ATAC_%A_%a.err
#SBATCH --partition=pibu_el8
#SBATCH --array=0-3
set -euo pipefail
 
echo $(date)
source "$SLURM_SUBMIT_DIR/00_config.sh"

/data/users/vmuller/miniforge3/envs/scReadSim/bin/Rscript --version

activate_conda
export RETICULATE_PYTHON="$CONDA_PREFIX/bin/python"


MANIFEST="/data/users/vmuller/0_master_thesis/data/SEAD_Dataset/models/ATAC/manifest.csv"
MODEL_DIR=$(awk -v line="$((SLURM_ARRAY_TASK_ID + 1))" -F',' 'NR==line {print $2}' "$MANIFEST")
GROUPNAME=$(awk -v line="$((SLURM_ARRAY_TASK_ID + 1))" -F',' 'NR==line {print $1}' "$MANIFEST")
SIM_OUT="/data/users/vmuller/0_master_thesis/data/SEAD_Dataset/patient_subsets/simulated_data/ATAC/$GROUPNAME"
 
mkdir -p "$SIM_OUT"
 
Rscript "$WORKDIR/src/Synthetic_data_generation/simulate_ATAC.R" \
  --modeldir  "$MODEL_DIR" \
  --outdir    "$SIM_OUT" \
  --donor_id  "$GROUPNAME"\
  --n_rep 50
 
echo "Finished"
echo $(date)

