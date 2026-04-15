#!/bin/bash
#SBATCH --job-name=run_mapmycells
#SBATCH --time=01:00:00
#SBATCH --mem=6G
#SBATCH --cpus-per-task=4
#SBATCH --output=/data/users/vmuller/0_master_thesis/log/run_mapmycells_%A_%a.out
#SBATCH --error=/data/users/vmuller/0_master_thesis/log/run_mapmycells_%A_%a.err
#SBATCH --partition=pshort_el8
#SBATCH --array=0-49%5

set -euo pipefail

# Load configuration and environment variables
source "$SLURM_SUBMIT_DIR/00_config.sh"

#Set up the conda environment for mapmycells
unset PYTHONPATH
export MAMBA_ROOT_PREFIX=$HOME/.conda
source ~/miniforge3/etc/profile.d/conda.sh
source ~/miniforge3/etc/profile.d/mamba.sh
conda activate cell_type_mapper

#Define the group and input directory for the simulated files
GROUP="Not AD"
SIMULATED_FILES_DIR="/data/users/vmuller/0_master_thesis/data/SEAD_Dataset/patient_subsets/simulated_data/RNA"
LIST_FILES=("$SIMULATED_FILES_DIR/$GROUP"/*.h5ad)

# Check if files exist and if the array index is within bounds
if [[ ! -e "${LIST_FILES[0]}" ]]; then
    echo "ERROR: No .h5ad files found in $SIMULATED_FILES_DIR/$GROUP/"
    exit 1
fi

echo "Total files found: ${#LIST_FILES[@]}"

if [[ $SLURM_ARRAY_TASK_ID -ge ${#LIST_FILES[@]} ]]; then
    echo "ERROR: Array index $SLURM_ARRAY_TASK_ID out of bounds (${#LIST_FILES[@]} files)"
    exit 1
fi

# Retrieve the file corresponding to the current array index and extract the base name and extension
file="${LIST_FILES[$SLURM_ARRAY_TASK_ID]}"
BASENAME=$(basename "$file")
FILENAME="${BASENAME%.*}"
EXTENSION="${BASENAME##*.}"

# Define paths for precomputed stats, temporary files, output results, and graphs
PATH_PRECOMPUTED_STATS="$WORKDIR/data/mapmycells/precomputed_stats/precomputed_stats.20231120.sea_ad.MTG.h5"
TMP="$WORKDIR/data/SEAD_Dataset/TEMP/mapmycells/$GROUP"
OUTPUTDIR="$WORKDIR/data/SEAD_Dataset/patient_subsets/mapmycells_results/$GROUP"
GRAPH_OUTPUTDIR="$WORKDIR/data/SEAD_Dataset/patient_subsets/mapmycells_results/$GROUP/graphs"

# Create necessary directories if they don't exist
mkdir -p "$TMP" "$OUTPUTDIR" "$GRAPH_OUTPUTDIR"

echo "Running mapmycells for group: $GROUP"
echo "Processing file: $file"

# Step 1: Convert to Ensembl ID format (skip if already done)
if [[ -f "$TMP/${FILENAME}_ens.$EXTENSION" ]]; then
    echo "Already converted, skipping: $FILENAME"
else
    python "$PATH_SCRIPTS/mapmycell_convert_genename_to_ENS.py" \
        "$file" \
        "$TMP"
fi

# Step 2: Run mapmycells (skip if already done)
if [[ -f "$OUTPUTDIR/$FILENAME.csv" ]]; then
    echo "Already mapped, skipping: $FILENAME"
else
    cd "$OUTPUTDIR"
    python -m cell_type_mapper.cli.map_to_on_the_fly_markers \
        --query_path "$TMP/${FILENAME}_ens.$EXTENSION" \
        --extended_result_path "$FILENAME.json" \
        --csv_result_path "$FILENAME.csv" \
        --n_processors 4 \
        --cloud_safe False \
        --precomputed_stats.path "$PATH_PRECOMPUTED_STATS" \
        --type_assignment.normalization raw \
        --query_markers.n_per_utility 15 \
        --reference_markers.log2_fold_min_th 0.5
fi

# Step 3: Generate graphs and summary statistics
python "$PATH_SCRIPTS/analysis_mapmycells.py" \
    "$FILENAME" \
    "$OUTPUTDIR/$FILENAME.csv" \
    "$GRAPH_OUTPUTDIR"

echo "Done: $FILENAME"