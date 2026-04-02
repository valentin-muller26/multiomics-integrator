#!/bin/bash
#SBATCH --job-name=download_mapmycells
#SBATCH --time=01:00:00
#SBATCH --mem=6G
#SBATCH --cpus-per-task=4
#SBATCH --output=/data/users/vmuller/0_master_thesis/log/download_mapmycells_%J.out
#SBATCH --error=/data/users/vmuller/0_master_thesis/log/download_mapmycells_%J.err
#SBATCH --partition=pshort_el8

set -euo pipefail

source "$SLURM_SUBMIT_DIR/00_config.sh"

#Installation of the map_my_cells software
unset PYTHONPATH
export MAMBA_ROOT_PREFIX=$HOME/.conda
source ~/miniforge3/etc/profile.d/conda.sh
source ~/miniforge3/etc/profile.d/mamba.sh
conda activate cell_type_mapper
pip install mygene
pip install scikit-learn
# Install cell_type_mapper if not already present
if ! python -c "import cell_type_mapper" 2>/dev/null; then
    cd "$WORKDIR/data/mapmycells"
    if [ ! -d "cell_type_mapper" ]; then
        git clone https://github.com/AllenInstitute/cell_type_mapper.git
    fi
    pip install -e cell_type_mapper/
fi



