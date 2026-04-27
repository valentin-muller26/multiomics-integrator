#!/bin/bash
#SBATCH --job-name=download_peaks
#SBATCH --time=1:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=/data/users/vmuller/0_master_thesis/log/download_peaks_%J.out
#SBATCH --error=/data/users/vmuller/0_master_thesis/log/download_peaks_%J.err
#SBATCH --partition=pshort_el8

set -euo pipefail
echo $(date)
source "$SLURM_SUBMIT_DIR/00_config.sh"

# Download peaks
PEAKS_DIR="/data/users/lfalquet/PYMantel/ATAC/VMuller/"
OUTPUT_DIR="$PATH_DATA/peaks_summits"
MERGE_TOOLS_DIR="$WORKDIR/data/ATAC_IterativeOverlapPeakMerging"

mkdir -p "$MERGE_TOOLS_DIR"
mkdir -p "$OUTPUT_DIR"

find "$PEAKS_DIR" -name "*summits.bed" -exec cp {} "$OUTPUT_DIR" \;


cd "$MERGE_TOOLS_DIR"
if [ ! -d "ATAC_IterativeOverlapPeakMerging" ]; then
    git clone https://github.com/corceslab/ATAC_IterativeOverlapPeakMerging.git
fi
