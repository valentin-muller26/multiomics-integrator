#!/bin/bash
#SBATCH --job-name=mapmycells_all_ADNC_group
#SBATCH --time=04:00:00
#SBATCH --mem=6G
#SBATCH --cpus-per-task=4
#SBATCH --output=/data/users/vmuller/0_master_thesis/log/mapmycells_all_ADNC_group_%J.out
#SBATCH --error=/data/users/vmuller/0_master_thesis/log/mapmycells_all_ADNC_group_%J.err
#SBATCH --partition=pibu_el8

set -euo pipefail

source "$SLURM_SUBMIT_DIR/00_config.sh"

unset PYTHONPATH
export MAMBA_ROOT_PREFIX="$HOME/miniforge3"
source ~/miniforge3/etc/profile.d/conda.sh
source ~/miniforge3/etc/profile.d/mamba.sh
conda activate cell_type_mapper

ADNC_STATUS=("Not AD" "Low" "Intermediate" "High")

TMP="$WORKDIR/data/SEAD_Dataset/TEMP/mapmycells"
mkdir -p $TMP

JOB1=$(sbatch --parsable 05b_run_mapmycell.sh "${ADNC_STATUS[0]}")
echo "Running Job for ${ADNC_STATUS[0]} group (Job ID: $JOB1)"

JOB2=$(sbatch --parsable --dependency=afterok:$JOB1 05b_run_mapmycell.sh "${ADNC_STATUS[1]}")
echo "Running Job for ${ADNC_STATUS[1]} group (Job ID: $JOB2)"

JOB3=$(sbatch --parsable --dependency=afterok:$JOB2 05b_run_mapmycell.sh "${ADNC_STATUS[2]}")
echo "Running Job for ${ADNC_STATUS[2]} group (Job ID: $JOB3)"

JOB4=$(sbatch --parsable --dependency=afterok:$JOB3 05b_run_mapmycell.sh "${ADNC_STATUS[3]}")
echo "Running Job for ${ADNC_STATUS[3]} group (Job ID: $JOB4)"

sbatch --dependency=afterok:$JOB4 --wrap="rm -rf $TMP" \
       --job-name=cleanup_mapmycells \
       --time=00:05:00 --mem=1G \
       --output=$WORKDIR/log/cleanup_%J.out \
       --partition=pibu_el8