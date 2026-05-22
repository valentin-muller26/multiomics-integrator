# MOSim pipeline

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

This step uses the scripts :
1. `02_simulation_mosim.sh` used for slurm to run the job and give the parameter to the R script `02_run_mosim.R`.
2. `02_run_mosim.R` :  script containing the pipeline step for the simulation of MOSim.
3. `mosim_functions.R` : script containing the different function used in the MOSim pipeline

**Input required**
The R script required the following input 
- donor_ID : prefix of the donor in the format ("DA16)
- inputdir : path to the folder containing the RNA-seq and ATAC-seq data extracted in the step 1
- outdir : path to the output directory where the results of the pipeline will be saved

**Output generated**

The MOSim simulation generates the following outputs in the `simulated_data/` folder:

- **`real/`**: real feature counts for the RNA-seq and ATAC-seq modalities (one file per donor and per modality), saved after the preprocessing step performed at the beginning of the simulation.
- **`merged/`**: simulated feature counts. For each donor, two files are produced:
  - one for RNA-seq, containing all simulated replicates,
  - one for ATAC-seq, containing all simulated replicates.

## Step 3 : Validation of the simulated data


## (optional) Step 4 : Sanity check to compare the distribution of two different donor real data

