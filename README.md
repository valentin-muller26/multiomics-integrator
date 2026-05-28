# MultiOmic analysis
This repository consists of the code and configuration for three pipelines developed as part of a Master's thesis in Bioinformatics and Computational Biology (Universities of Fribourg & Bern) on the multi-omics integration of bulk RNA-seq and ATAC-seq data from healthy donors in the context of atopic dermatitis.

## Repository overview

The repository is divided into three independent pipelines: 

1. **Single-cell dataset simulation using scDesign3** : In this pipeline, the `SEA-AD MTG Multiome` dataset was used as a proof-of-concept dataset to evaluate the feasibility of simulating paired single-cell RNA-seq and ATAC-seq data with `scDesign3` and consists of the following steps : downloading, preprocessing and quality control of the SEA-AD MTG Multiome dataset, simulation with `scDesign3`, and validation of the simulated counts. More details about the pipeline can be found in the  [pipeline README](src/scDesign_pipeline/README.md).

2. **Simulation of the bulk RNA-seq and ATAC-seq data from healthy donors dataset using MOSim** : Preprocessing of the bulk RNA-seq and ATAC-seq from healthy donors in the context of atopic dermatitis, simulation with `MOSim`, and validation with `countsimQC`. More details about the pipeline can be found in the [pipeline README](src/Mosim_pipeline/README.md).


3. **Multi-omics analysis of the bulk RNA-seq and ATAC-seq data from healthy donors dataset** :  Preprocessing of the bulk RNA-seq and ATAC-seq from healthy donors, training and downstream analysis of  unsupervised `MOFA2` models, and training of the supervised model `DIABLO` (mixOmics).  More details about the pipeline can be found in the [pipeline README](src/MultiOmic_analysis/README.md).

## Structure of the repository 

```
├── data            # Directory containing the datasets and the result of the different pipeline
├── envs            # Conda/mamba environment files (.yml)
├── log             # SLURM job logs, consisting of .err files containing errors and .out files containing the output of the jobs
├── README.md
└── src             # Directory containing the code of the three pipelines
    ├── scDesign_pipeline/
    ├── Mosim_pipeline/
    └── MultiOmic_analysis/
```


## Installation
For reproducibility, the exact versions of all R and Python dependencies are listed in separate `.yml` files (one per pipeline, with the exception of the `scDesign3` pipeline, which also requires the `cell_type_mapper` environment for the validation step).

Example of setting up an environment :
​

**Clone the repository**
```bash
git clone <repo-url>
cd <repo-name>
```

**Create the environment**

```bash
mamba env create -f envs/<environment>.yml
```

**Activate the environment**

```bash
mamba activate <environment>
```

### Requirements
- Linux operating system
- An HPC cluster managed by the SLURM workload manager
- `conda` or `mamba` for environment management (mamba recommended for
  solver speed)
- A minimum of 300 GB of memory is required for the simulation and
  validation steps of the `scDesign3` pipeline

## Running the pipelines using SLURM
All three pipelines were created to run on an HPC cluster using SLURM.
Each pipeline consists of a sequence of jobs. Each job is a bash script (`.sh`) that runs the corresponding Python or R script with the appropriate input parameters, and must be submitted in the order indicated by its numerical prefix.

Example of submitting job :

```bash

cd src/<pipeline>
sbatch 01_<name>.sh
sbatch 02_<name>.sh
...
```

## Author
Valentin Müller
