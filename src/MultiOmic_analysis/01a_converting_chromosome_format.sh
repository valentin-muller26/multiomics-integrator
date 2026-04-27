#!/bin/bash
#SBATCH --job-name=convert_chromosome_format
#SBATCH --time=1:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=/data/users/vmuller/0_master_thesis/log/convert_chromosome_format_%J.out
#SBATCH --error=/data/users/vmuller/0_master_thesis/log/convert_chromosome_format_%J.err
#SBATCH --partition=pshort_el8

set -euo pipefail
echo "$(date)"
source "$SLURM_SUBMIT_DIR/00_config.sh"
activate_conda

SUMMITS_DIR="$PATH_DATA/peaks_summits"
OUTPUT_DIR="$PATH_DATA/peaks_summits_converted"
mkdir -p "$OUTPUT_DIR"

for FILE in "$SUMMITS_DIR"/*_summits.bed; do
    OUTPUT_FILE="$OUTPUT_DIR/$(basename "$FILE")"
    awk 'BEGIN{OFS="\t"} {
        if ($1 == "MT") $1 = "chrM";
        else if ($1 !~ /^chr/) $1 = "chr" $1;
        print
    }' "$FILE" > "$OUTPUT_FILE"
done