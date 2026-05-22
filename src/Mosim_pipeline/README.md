# Mosim pipeline

## Pipeline overview and goals



## Step 0 : environment configuration
This step consists of two scripts:

1. **`00_config.sh`**: centralizes the main variables (paths to data,logs, working directory...) and the activation step of the conda/mamba environment.
2. **`00_setup_virtual_env.sh`**: creates or updata the conda/mamba virtual environment, installs the latest versions of all dependencies, and exports the list of dependencies to a `.yml` file in the `envs/` folder.

## Step 1 : Data extraction 
This step uses the script `01_data_extraction.sh`, which performs the following operations:

1. Copy the folder containing the quality control and feature count analyses produced by Prof. Falquet into the `data/` folder.
2. Extract the two RNA-seq feature count files from the copied folder and place them in `data_interleukine/RNA/`.
3. Extract the consensus feature counts and per-donor annotations for ATAC-seq from the copied folder and place them in `data_interleukine/ATAC/`.
4. Extract the consensus feature count for the combined ATAC-seq of donors HC16–HC19 and place it in `data_interleukine/ATAC_consensus/`.

## Step 2 : Simulation using MOsim

## Step 3 : Validation of the simulated data


## (optional) Step 4 : Sanity check to compare the distribution of two different donor real data

