#!/bin/bash
#SBATCH --job-name=setup_mosim_env
#SBATCH --time=1:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=/data/users/vmuller/0_master_thesis/log/mosim_setup_%J.out
#SBATCH --error=/data/users/vmuller/0_master_thesis/log/mosim_setup_%J.err
#SBATCH --partition=pshort_el8

set -euo pipefail
echo $(date)
source "$SLURM_SUBMIT_DIR/00_config.sh"

# Setup conda/mamba
unset PYTHONPATH
export MAMBA_ROOT_PREFIX=$HOME/.conda
source ~/miniforge3/etc/profile.d/conda.sh
source ~/miniforge3/etc/profile.d/mamba.sh

# 1. Create the conda environment if it does not exist
if ! conda env list | grep -q "^mosim "; then
    mamba create -n mosim -c conda-forge -c bioconda -y \
        r-base=4.3 \
        r-dplyr \
        r-argparse \
        bioconductor-mosim \
        bioconductor-deseq2 \
        bioconductor-countsimqc \
        python=3.10 pip
fi

# 2. Activate the environment
conda activate mosim


# 3. Export environment to yml
ENV_YML="$WORKDIR/envs/mosim.yml"
mkdir -p "$WORKDIR/envs"
conda env export > "$ENV_YML"
echo "Environment saved to ${ENV_YML}"

echo "Done. Activate with: conda activate mosim"