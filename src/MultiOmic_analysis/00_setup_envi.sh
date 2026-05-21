#!/bin/bash
#SBATCH --job-name=setup_multiomic_analysis_env
#SBATCH --time=2:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=/data/users/vmuller/0_master_thesis/log/setup_multiomic_analysis_env_%J.out
#SBATCH --error=/data/users/vmuller/0_master_thesis/log/setup_multiomic_analysis_env_%J.err
#SBATCH --partition=pshort_el8

set -euo pipefail
echo "$(date)"
source "$SLURM_SUBMIT_DIR/00_config.sh"

# Setup conda/mamba
unset PYTHONPATH
export MAMBA_ROOT_PREFIX=$HOME/.conda
source ~/miniforge3/etc/profile.d/conda.sh
source ~/miniforge3/etc/profile.d/mamba.sh

# List of packages to install
PACKAGES=(
    r-base=4.3
    r-optparse
    r-argparse
    r-dplyr
    r-purrr
    r-stringr
    r-reshape2
    r-tidyr
    r-yaml
    r-readr
    r-data.table
    r-ggplot2
    r-ggrepel
    r-pheatmap
    r-rcolorbrewer
    r-cowplot
    r-psych
    r-mass
    r-matrixstats
    r-openxlsx
    r-ggpubr
    bioconductor-rtracklayer
    bioconductor-biostrings
    bioconductor-summarizedexperiment
    bioconductor-genomeinfodb
    bioconductor-genomicranges
    bioconductor-bsgenome
    bioconductor-bsgenome.hsapiens.ucsc.hg38
    bioconductor-deseq2
    bioconductor-edger
    bioconductor-limma
    bioconductor-mofa2
    bioconductor-mixomics
    bioconductor-mofadata
    bioconductor-org.hs.eg.db
    bioconductor-annotationdbi
    bioconductor-biomart
    bioconductor-clusterprofiler
    bioconductor-txdb.hsapiens.ucsc.hg38.knowngene
    r-msigdbr 
    bioconductor-fgsea
    bioconductor-reactomepa
    bioconductor-chipseeker
    python=3.10
    pip
    numpy
    scipy
    pandas
    h5py
    scikit-learn
    mofapy2
)

# 1. Create environment if absent, otherwise install/update missing packages
if conda env list | grep -q "^multiomic_analysis "; then
    echo "Environment 'multiomic_analysis' exists — installing/updating packages..."
    mamba install -n multiomic_analysis -c conda-forge -c bioconda -y "${PACKAGES[@]}"
else
    echo "Creating environment 'multiomic_analysis'..."
    mamba create -n multiomic_analysis -c conda-forge -c bioconda -y "${PACKAGES[@]}"
fi

# 2. Activate environment
conda activate multiomic_analysis



# 3. Export environment to yml
ENV_YML="$WORKDIR/envs/multiomic_analysis.yml"
mkdir -p "$WORKDIR/envs"
conda env export > "$ENV_YML"
echo "Environment saved to ${ENV_YML}"

echo "Done. Activate with: conda activate multiomic_analysis"