# Multi-omics analysis pipeline


This folder contains two versions of the MOFA2 pipeline and one version of the training step of the tool Diablo : a local version, intended to be run on a personal computer, and a cluster version, intended to be run on an HPC cluster via SLURM. Both versions implement the same MOFA2 workflow; the cluster version adds SLURM submission scripts and environment management.

The MOFA2 pipeline consists of two steps:
1. Preprocessed the raw RNA-seq and ATAC-seq datasets and Train the model
2. Downstream analysis : generates the visualization plots, performs Gene Set Enrichment Analysis (GSEA), and links ATAC-seq peaks to their nearest genes using CHIPSeeker.

# Local version

The local version consists of two helper files and two scripts corresponding to the two steps of the pipeline:

- `utils.R` : contains functions for preprocessing, normalization, batch correction, and feature selection via Median Absolute Deviation (MAD), as well as functions for generating diagnostic PCA plots. These functions are shared between the MOFA2 and DIABLO pipelines.
- `MOFA_function.R` : contains functions specific to the MOFA2 pipeline (training, GSEA, etc.).
- `MOFA_train.R` : performs preprocessing, normalization, feature selection, MOFA2 model training, and selection of the best model.
- `MOFA_downstream_analysis.R` : runs the downstream analysis step.

## Requirements

To run this pipeline, download the MOFA_local/ folder along with the RNA-seq dataset and the ATAC-seq consensus calls. Rename them to RNAseq.txt and ATACseq.txt respectively, and place them in a data/ folder inside the folder containing the scripts.

# Cluster version 

In addition to the adapted R scripts from the local version, the cluster version provides bash scripts to run the different steps of the MOFA pipeline, along with the following configuration and setup scripts:
- `00_config.sh` : centralizes the main variables (paths to data,logs, working directory...) and the activation step of the conda/mamba environment.
 - `00_setup_envi.sh` : creates or update the conda/mamba virtual environment, installs the latest versions of all dependencies, and exports the list of dependencies to a .yml file in the envs/ folder.
- `01_data_extraction.sh` :  which performs the following operations:
  - Copy the folder containing the quality control and feature count analyses produced by Prof. Falquet into the data/ folder.
  - Extract the two RNA-seq feature count files from the copied folder and place them in data_interleukine/RNA/.
  - Extract the consensus feature counts and per-donor annotations for ATAC-seq from the copied folder and place them in data_interleukine/ATAC/.
  - Extract the consensus feature count for the combined ATAC-seq of donors HC16–HC19 and place it in data_interleukine/ATAC_consensus/.
