#!/bin/bash
#SBATCH --job-name=setup_mapmycells_env
#SBATCH --time=01:00:00
#SBATCH --mem=6G
#SBATCH --cpus-per-task=4
#SBATCH --output=/data/users/vmuller/0_master_thesis/log/setup_mapmycells_env_%J.out
#SBATCH --error=/data/users/vmuller/0_master_thesis/log/setup_mapmycells_env_%J.err
#SBATCH --partition=pshort_el8

set -euo pipefail

source "$SLURM_SUBMIT_DIR/00_config.sh"

# Setup conda/mamba
unset PYTHONPATH
export MAMBA_ROOT_PREFIX="$HOME/miniforge3"
source ~/miniforge3/etc/profile.d/conda.sh
source ~/miniforge3/etc/profile.d/mamba.sh

# 1. Create and activate cell_type_mapper environment
if ! conda env list | awk '{print $1}' | grep -qx "cell_type_mapper"; then
    echo "Creating cell_type_mapper environment..."
    mamba create -n cell_type_mapper python=3.10 -y
fi

conda activate cell_type_mapper

# Install pip deps only if missing
if ! python -c "import mygene, sklearn" 2>/dev/null; then
    pip install mygene scikit-learn
fi

# 2. Install cell_type_mapper package if not already installed
if ! python -c "import cell_type_mapper" 2>/dev/null; then
    mkdir -p "$WORKDIR/data/mapmycells"
    cd "$WORKDIR/data/mapmycells"
    if [ ! -d "cell_type_mapper" ]; then
        git clone https://github.com/AllenInstitute/cell_type_mapper.git
    fi
    pip install -e cell_type_mapper/
fi

# 3. Download SEA-AD MTG precomputed stats file
PRECOMPUTED_DIR="$WORKDIR/data/mapmycells/precomputed_stats"
mkdir -p "$PRECOMPUTED_DIR"
if [ ! -f "$PRECOMPUTED_STATS_FILE" ]; then
    echo "Downloading precomputed stats..."
    wget -c -O "$PRECOMPUTED_STATS_FILE" "$PRECOMPUTED_STATS_URL"
else
    echo "Precomputed stats already exist, skipping."
fi

# 4. Export environment to yml
ENV_YML="$WORKDIR/envs/cell_type_mapper.yml"
mkdir -p "$WORKDIR/envs"
conda env export > "$ENV_YML"
echo "Environment saved to ${ENV_YML}"