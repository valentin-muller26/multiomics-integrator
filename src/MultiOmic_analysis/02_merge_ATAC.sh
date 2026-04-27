#!/bin/bash
#SBATCH --job-name=merge_peaks
#SBATCH --time=2:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=/data/users/vmuller/0_master_thesis/log/merge_peaks_%J.out
#SBATCH --error=/data/users/vmuller/0_master_thesis/log/merge_peaks_%J.err
#SBATCH --partition=pshort_el8

# 
set -euo pipefail
echo $(date)
source "$SLURM_SUBMIT_DIR/00_config.sh"
activate_conda


# Define and create necessary directories
MERGE_TOOLS_DIR="$WORKDIR/data/ATAC_IterativeOverlapPeakMerging/ATAC_IterativeOverlapPeakMerging"
SUMMITS_DIR="$PATH_DATA/peaks_summits_converted/"
#SUMMITS_DIR="$PATH_DATA/test_peak/"
OUTPUT_DIR="$PATH_DATA/merged_peaks/"
METADATA_DIR="$SUMMITS_DIR/metadata/"
BLACKLIST_DIR="$PATH_DATA/blacklist"
BLACKLIST_FILE="$BLACKLIST_DIR/hg38.blacklist.bed"
mkdir -p "$BLACKLIST_DIR"
mkdir -p "$METADATA_DIR"
mkdir -p "$OUTPUT_DIR"


# Download blacklist file for hg38 if not already present
if [ ! -f "$BLACKLIST_FILE" ]; then
    echo "Downloading hg38 blacklist..."
    cd "$BLACKLIST_DIR"
    wget -q http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/hg38.blacklist.bed.gz
    gunzip hg38.blacklist.bed.gz
    cd -
fi


# Remove old metadata file if exists
if [ -f "$METADATA_DIR/metadata.txt" ]; then
    rm "$METADATA_DIR/metadata.txt"
fi

echo -e "Sample\tGroup" > "$METADATA_DIR/metadata.txt"
# Create metadata file for the merging tool
ls "$SUMMITS_DIR"/*_summits.bed | while read -r file
do
    filename=$(basename "$file")
    sample_name="${filename%_summits.bed}"
    group=$(echo "$sample_name" | cut -d'_' -f2)
    echo -e "${sample_name}\t${group}" >> "$METADATA_DIR/metadata.txt"
done


# Run the merging tool
Rscript "$MERGE_TOOLS_DIR/createIterativeOverlapPeakSet.R" \
    --metadata "$METADATA_DIR/metadata.txt" \
    --macs2dir "$SUMMITS_DIR" \
    --outdir "$OUTPUT_DIR" \
    --blacklist "$BLACKLIST_DIR/hg38.blacklist.bed" \
    --suffix _summits.bed \
    --genome hg38 \
    --spm 5 \
    --rule "(n+1)/2" \
    --extend 250
       