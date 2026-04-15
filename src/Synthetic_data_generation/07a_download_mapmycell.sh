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

# Setup conda/mamba
unset PYTHONPATH
export MAMBA_ROOT_PREFIX=$HOME/.conda
source ~/miniforge3/etc/profile.d/conda.sh
source ~/miniforge3/etc/profile.d/mamba.sh

# 1. Create environment if not already present
if ! conda env list | grep -q "^cell_type_mapper "; then
    echo "Creating cell_type_mapper environment..."
    mamba create -n cell_type_mapper python=3.10 -y
fi

# 2. Activate and install dependencies
conda activate cell_type_mapper

pip install mygene scikit-learn

if ! python -c "import cell_type_mapper" 2>/dev/null; then
    mkdir -p "$WORKDIR/data/mapmycells"
    cd "$WORKDIR/data/mapmycells"
    if [ ! -d "cell_type_mapper" ]; then
        git clone https://github.com/AllenInstitute/cell_type_mapper.git
    fi
    pip install -e cell_type_mapper/
fi

# 3. Export environment to yml
ENV_YML="$WORKDIR/envs/cell_type_mapper.yml"
mkdir -p "$WORKDIR/envs"
conda env export > "$ENV_YML"
echo "Environment saved to ${ENV_YML}"