# scDesign pipeline

## Pipeline overview and goals


## Step 0 : environment configuration
This step consists of two scripts:

1. **`00_config.sh`**: centralizes the main variables (paths to data,logs, working directory, URL of the AWS S3 repository, ...) and the activation step of the conda/mamba environment.
2. **`00_setup_scDesign_env.sh`**: creates or update the conda/mamba virtual environment, installs the latest versions of all dependencies, and exports the list of dependencies to a `.yml` file in the `envs/` folder.

# Step 1 Downloading the dataset
This step consists of the script `01_download_SEAD_Dataset.sh`, which downloads the final 2024-02-13 release of the single-cell RNA-seq and single-cell ATAC-seq data from the SEA-AD consortium's official AWS S3 repository.

# Step 2 : Dataset preprocessing 
