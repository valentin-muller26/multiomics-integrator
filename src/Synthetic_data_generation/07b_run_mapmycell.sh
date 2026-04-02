#!/bin/bash
#SBATCH --job-name=run_mapmycells
#SBATCH --time=01:00:00
#SBATCH --mem=6G
#SBATCH --cpus-per-task=4
#SBATCH --output=/data/users/vmuller/0_master_thesis/log/run_mapmycells_%J.out
#SBATCH --error=/data/users/vmuller/0_master_thesis/log/run_mapmycells_%J.err
#SBATCH --partition=pshort_el8
#SBATCH --array=0-49

set -euo pipefail

source "$SLURM_SUBMIT_DIR/00_config.sh"

#Installation of the map_my_cells software
unset PYTHONPATH
export MAMBA_ROOT_PREFIX=$HOME/.conda
source ~/miniforge3/etc/profile.d/conda.sh
source ~/miniforge3/etc/profile.d/mamba.sh
conda activate cell_type_mapper


#Simulated files to use for mapmycells
GROUP="NOT_AD"
SIMULATED_FILES_DIR="/data/users/vmuller/0_master_thesis/data/SEAD_Dataset/patient_subsets/simulated_data"
LIST_FILES=("$SIMULATED_FILES_DIR"/$GROUP/*.h5ad)
file="${LIST_FILES[$SLURM_ARRAY_TASK_ID]}"


BASENAME=$(basename "$file") #Retrieve the file name with extension (ex: sample_001.h5ad)        
FILENAME="${BASENAME%.*}"    #Remove the extension to get the file name without extension (ex: sample_001)
EXTENSION="${BASENAME##*.}"  #Retrieve the extension (ex: h5ad)

PATH_PRECOMPUTED_STATS="$WORKDIR/data/mapmycells/precomputed_stats/precomputed_stats.20231120.sea_ad.MTG.h5"
TMP="$WORKDIR/data/SEAD_Dataset/TEMP/mapmycells/$GROUP"
OUTPUTDIR="$WORKDIR/data/SEAD_Dataset/patient_subsets/mapmycells_results/$GROUP"
GRAPH_OUTPUTDIR="$WORKDIR/data/SEAD_Dataset/patient_subsets/mapmycells_results/$GROUP/graphs"

mkdir -p "$OUTPUTDIR"
mkdir -p "$TMP"
mkdir -p "$GRAPH_OUTPUTDIR"

echo "Running mapmycells for group $GROUP"
echo "Running mapmycells for file: $file"

# Convert the simulated files to ENS format
python "$WORKDIR/src/Synthetic_data_generation/mapmycell_convert_genename_to_ENS.py" \
   $file \
   $TMP

#Run mapmycells
cd $OUTPUTDIR
python -m cell_type_mapper.cli.map_to_on_the_fly_markers \
    --query_path "$TMP/${FILENAME}_ens.$EXTENSION" \
    --extended_result_path "$FILENAME.json" \
    --csv_result_path "$FILENAME.csv" \
    --n_processors 1 \
    --cloud_safe False \
    --precomputed_stats.path "$PATH_PRECOMPUTED_STATS" \
    --type_assignment.normalization raw \
    --query_markers.n_per_utility 15 \
    --reference_markers.log2_fold_min_th 0.5 


#Create the graphs for the results of mapmycells and the summary statistics
GRAPHDIR="$GRAPH_OUTPUTDIR"
python "$WORKDIR/src/Synthetic_data_generation/analysis_mapmycells.py" \
    "$FILENAME" \
    "$OUTPUTDIR/$FILENAME.csv" \
    "$GRAPHDIR"