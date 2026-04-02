#!/bin/bash
#SBATCH --job-name=extracting_data
#SBATCH --time=1:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --output=/data/users/vmuller/0_master_thesis/log/extracting_data_%J.out
#SBATCH --error=/data/users/vmuller/0_master_thesis/log/extracting_data_%J.err
#SBATCH --partition=pibu_el8

set -euo pipefail

#Define the variables for the directories and required files
OUTPUTDIR="/data/users/vmuller/0_master_thesis/data/data_interleukines"
INPUTDIR="/data/users/vmuller/0_master_thesis/data/VMuller"

#Create the output directory if it does not exist
mkdir -p "$OUTPUTDIR/RNA"
mkdir -p "$OUTPUTDIR/ATAC"

#Retrive the RNA-seq data
cp $INPUTDIR/*geneCounts.txt $OUTPUTDIR/RNA/


mapfile -t list_ATACSEQ_donors < <(find "$INPUTDIR/" -type d -maxdepth 1 -mindepth 1)

for donor in "${list_ATACSEQ_donors[@]}"; do
    if [[ -d "$donor/consensus" ]]; then
        cp "$donor/consensus/consensus_peaks.mLb.clN.annotatePeaks.txt" "$OUTPUTDIR/ATAC/${donor##*/}_consensus_peaks_annotatePeaks.txt"
        cp "$donor/consensus/consensus_peaks.mLb.clN.featureCounts.txt" "$OUTPUTDIR/ATAC/${donor##*/}_consensus_peaks_featureCounts.txt"
    fi
done

echo "Data extraction completed successfully."