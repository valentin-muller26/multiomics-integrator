#!/bin/bash
#SBATCH --job-name=setup_scdesign_env
#SBATCH --time=2:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --output=/data/users/vmuller/0_master_thesis/log/setup_scdesign_env_%J.out
#SBATCH --error=/data/users/vmuller/0_master_thesis/log/setup_scdesign_env_%J.err
#SBATCH --partition=pshort_el8

set -euo pipefail
echo "Start: $(date)"
source "$SLURM_SUBMIT_DIR/00_config.sh"

# --- Conda/Mamba setup ---
unset PYTHONPATH
export MAMBA_ROOT_PREFIX=$HOME/.conda
source ~/miniforge3/etc/profile.d/conda.sh
source ~/miniforge3/etc/profile.d/mamba.sh

ENV_NAME="scdesign"

# List of packages to install
PACKAGES=(
    r-base=4.3
    r-devtools
    r-biocmanager
    r-optparse
    r-dplyr
    r-irlba
    r-matrix
    r-ggplot2
    r-tidyr
    r-reshape2
    r-umap
    r-tibble
    r-anndata
    r-cowplot
    r-pbmcapply
    r-rfast
    r-scales
    r-ggrastr
    r-argparse
    bioconductor-cellmixs
    bioconductor-zellkonverter
    bioconductor-singlecellexperiment
    bioconductor-sparsematrixstats
    bioconductor-matrixgenerics
    bioconductor-scuttle
    bioconductor-deseq2
    python=3.9
    pip
    numpy
    pandas
    matplotlib
    seaborn
    anndata
)

# --- Create or update environment ---
if mamba env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
    echo "Environment '${ENV_NAME}' exists — updating packages..."
    mamba install -n "${ENV_NAME}" -c conda-forge -c bioconda -y "${PACKAGES[@]}"
else
    echo "Creating environment '${ENV_NAME}'..."
    mamba create -n "${ENV_NAME}" -c conda-forge -c bioconda -y "${PACKAGES[@]}"
fi

# --- Activate ---
conda activate "${ENV_NAME}"

# --- Install scDesign3 from GitHub ---
echo "Installing scDesign3 from GitHub..."
Rscript -e '
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# Make sure Bioconductor deps resolve correctly
options(repos = BiocManager::repositories())
if (!requireNamespace("scDesign3", quietly = TRUE)) {
     devtools::install_github("SONGDONGYUAN1994/scDesign3", upgrade = "never", dependencies = TRUE)
   }
# Sanity check
library(scDesign3)
cat("scDesign3 version:", as.character(packageVersion("scDesign3")), "\n")
'

# --- Export environment ---
ENV_DIR="$WORKDIR/envs"
mkdir -p "$ENV_DIR"
conda env export                 > "$ENV_DIR/${ENV_NAME}.yml"
echo "Environment exported to ${ENV_DIR}/${ENV_NAME}.yml"

echo "Done: $(date)"
echo "Activate with: conda activate ${ENV_NAME}"