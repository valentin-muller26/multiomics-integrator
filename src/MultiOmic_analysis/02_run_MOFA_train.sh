#!/bin/bash
#SBATCH --job-name=train_MOFA2
#SBATCH --time=1:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=4
#SBATCH --output=/data/users/vmuller/0_master_thesis/log/train_MOFA2_%J.out
#SBATCH --error=/data/users/vmuller/0_master_thesis/log/train_MOFA2_%J.err
#SBATCH --partition=pibu_el8

set -euo pipefail

echo "$(date)"

# Load the configuration and activate the conda environment
source "$SLURM_SUBMIT_DIR/00_config.sh"
activate_conda

# Define the variables for the input and output
INPUT_RNA_FILE="$PATH_DATA/RNA/RNAseq_run3_AO1_HC14_HC19_geneCounts.txt"
INPUT_ATAC_FILE="$PATH_DATA/ATAC_consensus/consensus_peaks.mLb.clN.featureCounts.txt"

# Create the output directories if not present
mkdir -p "$OUTPUT_MODEL" "$OUTPUT_GRAPH"

# Run the training pipeline
Rscript "$PATH_SCRIPTS/02_MOFA_train.R" \
  --nb_factor       5 \
  --input_path_rna  "$INPUT_RNA_FILE" \
  --input_path_atac "$INPUT_ATAC_FILE" \
  --outdir_model    "$OUTPUT_MODEL" \
  --outdir_graph    "$OUTPUT_GRAPH"

echo "MOFA training finished"
echo "$(date)"