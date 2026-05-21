#!/bin/bash
#SBATCH --job-name=MOFA2_downstream_analyis
#SBATCH --time=1:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=4
#SBATCH --output=/data/users/vmuller/0_master_thesis/log/MOFA2_downstream_analyis_%J.out
#SBATCH --error=/data/users/vmuller/0_master_thesis/log/MOFA2_downstream_analyis_%J.err
#SBATCH --partition=pibu_el8

set -euo pipefail

echo "$(date)"

# Load the configuration and activate the conda environment
source "$SLURM_SUBMIT_DIR/00_config.sh"
activate_conda


# Create the output directories if not present
mkdir -p "$OUTPUT_MODEL" "$OUTPUT_GRAPH"

# Run the training pipeline
Rscript "$PATH_SCRIPTS/03_MOFA_downstream_analysis.R" \
  --outdir_model    "$OUTPUT_MODEL" \
  --outdir_graph    "$OUTPUT_GRAPH"

echo "MOFA downstream analysis finished"
echo "$(date)"