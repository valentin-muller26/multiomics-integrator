#!/bin/bash
#SBATCH --job-name=extracting_data
#SBATCH --time=1:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --output=/data/users/vmuller/0_master_thesis/log/extracting_data_%J.out
#SBATCH --error=/data/users/vmuller/0_master_thesis/log/extracting_data_%J.err
#SBATCH --partition=pibu_el8

echo "$(date)"
set -euo pipefail
source "$SLURM_SUBMIT_DIR/00_config.sh"


# Define directories
SOURCE_DIR="/data/users/lfalquet/PYMantel/ATAC/VMuller"
OUTPUTDIR="$PATH_DATA"
INPUTDIR="$WORKDIR/data/VMuller"
MULTI_DONOR_ATAC_CONSENSUS_FILE="$INPUTDIR/ALL16-19_narrowpeaks/merged_library/consensus/consensus_peaks.mLb.clN.featureCounts.txt"

# Copy the data from the source directory to the working directory
echo "Copying data from $SOURCE_DIR to $WORKDIR/data/"
cp -r "$SOURCE_DIR" "$WORKDIR/data/"

#Create the output directory if it does not exist
mkdir -p "$OUTPUTDIR/RNA"
mkdir -p "$OUTPUTDIR/ATAC"
mkdir -p "$OUTPUTDIR/ATAC_consensus"

#Retrive the RNA-seq data
echo "Copying RNA-seq data"
cp "$INPUTDIR"/*geneCounts.txt "$OUTPUTDIR/RNA/"

# Retrieve the ATAC-seq data consensus for merged donor HC 16-19
echo "Copying ATAC-seq consensus data for merged donors HC 16-19"
cp "$MULTI_DONOR_ATAC_CONSENSUS_FILE" "$OUTPUTDIR/ATAC_consensus/"

# List ATAC-seq donor directories and copy consensus peaks files
echo "Copying ATAC-seq consensus peaks files for each donor"
mapfile -t list_ATACSEQ_donors < <(find "$INPUTDIR" -mindepth 1 -maxdepth 1 -type d)
echo "Found ${#list_ATACSEQ_donors[@]} ATAC donor directories"

for donor in "${list_ATACSEQ_donors[@]}"; do
    donor_name="${donor##*/}"

    # Skip the multi-donor merged directory (handled separately above)
    if [[ "$donor_name" == "ALL16-19_narrowpeaks" ]]; then
        echo "  - $donor_name: skipped (multi-donor, handled separately)"
        continue
    fi

    if [[ -d "$donor/consensus" ]]; then
        echo "  - $donor_name: copying consensus files"
        cp "$donor/consensus/consensus_peaks.mLb.clN.annotatePeaks.txt" \
           "$OUTPUTDIR/ATAC/${donor_name}_consensus_peaks_annotatePeaks.txt"
        cp "$donor/consensus/consensus_peaks.mLb.clN.featureCounts.txt" \
           "$OUTPUTDIR/ATAC/${donor_name}_consensus_peaks_featureCounts.txt"
    else
        echo "  - $donor_name: no consensus directory, skipping"
    fi
done


echo "Data extraction completed successfully."
echo "$(date)"