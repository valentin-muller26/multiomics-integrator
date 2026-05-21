#' Generate synthetic count replicates per ADNC group using scDesign3 (ATAC-seq)
#'
#' @description
#' Loads the manifest (TSV listing the .rds path of every batch model fitted
#' for one ADNC group by the previous script) and iterates over the requested
#' number of replicates. Each replicate is simulated using a distinct seed.
#' For each replicate every batch is simulated independently
#' the resulting count matrices are merged on their common peaks,
#' and the merged result is exported as a single .h5ad file.
#'
#' The `--start_replicate_count` argument shifts both the replicate index
#' and the seed, allowing an interrupted run to be resumed without overwriting
#' existing replicates.
#'
#' @param modeldir Path to the directory containing model_manifest.tsv and the per-batch fitted model (.rds) files.
#' @param outdir Path to the output directory where the simulated replicates (.h5ad) will be written. Created if it does not exist.
#' @param ADNC_cat Name of the ADNC category being simulated. Used as the filename prefix and stored in obs metadata of every simulated cell.
#' @param n_rep Number of replicates to generate [default: 10].
#' @param start_replicate_count Starting index for replicate numbering and seeding. Set > 1 to resume an interrupted run [default: 1].
#' @param family Distribution family for the marginal model: nb, zinb, or zip [default: 'nb']. Must match the family used at fitting time.
#'
#' @return Writes one .h5ad file per replicate to outdir, named <ADNC_cat>_rep<NN>_seed<S>_scDesign3.h5ad.



library(Matrix)
library(zellkonverter)
library(SingleCellExperiment)
library(anndata)
library(argparse)
library(scDesign3)

#  Parse arguments
parser <- ArgumentParser(description = "scDesign3 ATAC-seq — generate replicates from saved models")
parser$add_argument("-m", "--modeldir", type = "character", required = TRUE,
                    help = "Directory containing saved .rds model files")
parser$add_argument("-o", "--outdir",   type = "character", required = TRUE,
                    help = "Output directory for simulated h5ad replicates")
parser$add_argument("--ADNC_cat",       type = "character", default = "donor",
                    help = "ADNC identifier added to obs metadata [default: 'donor']")
parser$add_argument("--n_rep",          type = "integer",   default = 10,
                    help = "Number of replicates to generate [default: 10]")
parser$add_argument("--start_replicate_count", type = "integer", default = 1,
                    help = "Starting count for replicate numbering used to resume simulation [default: 1]")
parser$add_argument("--family",         type = "character", default = "nb",
                    choices = c("nb", "zinb", "zip"),
                    help = "Distribution family for the marginal model: nb (negative binomial), zinb (zero-inflated NB), or zip (zero-inflated Poisson) [default: 'nb']. Must match the family used at fitting time.")
args <- parser$parse_args()

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Simulation function that takes a fitted model object, a random seed, and a batch ID to generate new counts and metadata
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
simulate_from_model <- function(fit_obj, seed, batch_id, family) {
    #' Function to simulate new counts for 1 batch
    #'
    #' @description
    #' This function extracts marginal parameters from the fitted model and runs
    #' the simulation with the given seed and distribution family, renames the
    #' cells in the format (cell type + unique ID + batch ID), and builds the
    #' metadata for the obs field of the merged data. Errors are caught and
    #' logged, returning NULL so the calling loop can skip the batch without
    #' aborting the replicate.
    #'
    #' @param fit_obj : object containing the fitted model of one batch loaded from .rds
    #' @param seed : seed for the simulation (different per replicate)
    #' @param batch_id : number of the current batch, used for logging and cell naming
    #' @param family : distribution family for the marginal model ("nb", "zinb", "zip")
    #'
    #' @return : list containing the simulated counts and the metadata, or NULL if the simulation failed

    set.seed(seed)

    # Extract the new covariate data used for simulation from the fitted model object
    # New covariate contains the metadata for the new simulated cells, including cell type and donor id as defined during model fitting
    new_cov <- fit_obj$dat$newCovariate

    # Print dimensions and cell type distribution of the new covariate for understanding the input to the simulation
    # Rows correspond to cells, columns correspond to covariates containing cell type and donor id
    cat("    new_cov dims       :", nrow(new_cov), "cells x", ncol(new_cov), "(cell type, donor id)\n")
    if ("cell_type" %in% colnames(new_cov)) {
        # Retrieve the distribution of cell types in the new covariate and print it for verification
        cell_type_list <- table(new_cov$cell_type) # Count the number of cells for each cell type in the new covariate
        cat("    cell_type distribution  :",
            paste(names(cell_type_list), cell_type_list, sep = "=", collapse = ", "), "\n")
    }

    # Step 1 — extract parameters (deterministic given model + covariate)
    cat("    Extracting parameters...\n")
    para <- tryCatch(
        extract_para(
            sce           = fit_obj$sce_ref,
            marginal_list = fit_obj$marginal_fit,
            n_cores       = 4,
            family_use    = family,
            new_covariate = new_cov,
            data          = fit_obj$dat$dat
        ),
        error = function(e) {
            cat("    ERROR in the parameter extraction:", conditionMessage(e), "\n")
            NULL
        }
    )
    if (is.null(para)) return(NULL) # If parameter extraction fails, return NULL to skip this batch

    # Print the number of peaks and cells for which parameters were extracted to understand the scale of the simulation
    cat(sprintf("    Fitted model: %d cells x %d peaks\n", nrow(para$mean_mat), ncol(para$mean_mat)))

    # Step 2 — simulate counts using the extracted parameters and the fitted copula model
    cat("    Simulating counts...\n")
    counts <- tryCatch(
        simu_new(
            sce               = fit_obj$sce_ref,
            mean_mat          = para$mean_mat,
            sigma_mat         = para$sigma_mat,
            zero_mat          = para$zero_mat,
            quantile_mat      = NULL,
            copula_list       = fit_obj$copula_fit$copula_list,
            n_cores           = 4,
            family_use        = family,
            input_data        = fit_obj$dat$dat,
            new_covariate     = new_cov,
            important_feature = fit_obj$copula_fit$important_feature,
            filtered_gene     = fit_obj$dat$filtered_gene
        ),
        error = function(e) {
            cat("    ERROR in the simulation:", conditionMessage(e), "\n")
            NULL
        }
    )

    # If simulation fails or returns 0 cells, return NULL to skip this batch
    if (is.null(counts)) return(NULL)
    if (ncol(counts) == 0) {
        cat("    WARNING: Simulation returned 0 cells — skipping.\n")
        return(NULL)
    }

    # Print the dimensions of the simulated counts matrix (should correspond to peaks x simulated cells)
    cat("    Simulated counts  :", nrow(counts), "peaks x", ncol(counts), "cells\n")

    # Rename the columns of the counts matrix to include cell type, a unique identifier, and batch information for downstream analysis and merging
    number_cell_simulated <- ncol(counts)
    cell_type_vector      <- as.character(new_cov$cell_type[seq_len(number_cell_simulated)])
    counter               <- ave(seq_along(cell_type_vector), cell_type_vector, FUN = seq_along) # Create a counter to have unique ID for each cell within the same cell type
    cell_name             <- paste0(cell_type_vector, "_", counter, "_b", batch_id) # Name of the cell : cell type + unique ID + batch ID
    colnames(counts)      <- cell_name

    # Build metadata to have the cell_type in the obs of the h5ad file
    meta           <- new_cov[seq_len(number_cell_simulated), , drop = FALSE]
    rownames(meta) <- cell_name

    list(counts = counts, meta = meta)
}
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Main pipeline
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Main variable for the simulation
N_REP      <- args$n_rep
SEED_START <- args$start_replicate_count

# Create output directory if it doesn't exist
dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)

# Record start time and print parameters of the simulation
start_time <- Sys.time()
cat("ADNC_cat  :", args$ADNC_cat, "\n")
cat("n_rep     :", N_REP, "\n")
cat("family    :", args$family, "\n\n")

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  Load model files from the specified directory, error if no manifest present
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------
manifest_path <- file.path(args$modeldir, "model_manifest.tsv")
if (!file.exists(manifest_path)) stop("model_manifest.tsv not found in: ", args$modeldir)

#Read the manifest : the manifest contain the path to the different rds of the fitted model for each batch of the ADNC group
manifest  <- read.table(manifest_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
rds_files <- manifest$rds_path
cat("Loaded manifest:", manifest_path, "\n")
cat("Found", length(rds_files), "model file(s).\n\n")
if (length(rds_files) == 0) stop("No model files found in: ", args$modeldir)

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  Generate replicates
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
for (rep_k in seq_len(N_REP)) {
    #Increment the seed to have different simulation for each replicate
    seed_k <- SEED_START + rep_k -1
    cat(sprintf("\n======  Replicate %d/%d  (seed=%d)  ======\n",
                rep_k, N_REP, seed_k))

    # Pre-allocate a list to store results for each batch within this replicate
    replicat_results <- vector("list", length(rds_files))

    # Iterate over each model batch
    for (batch in seq_along(rds_files)) {
        rds_path <- rds_files[[batch]]
        cat(sprintf("  [Batch %d/%d] Loading: %s\n", batch, length(rds_files),
                    basename(rds_path)))
        t0 <- Sys.time()

        # Load the object containing the values of the fitted model
        fit_model <- tryCatch(
            readRDS(rds_path),
            error = function(e) {
                cat("    ERROR loading model:", conditionMessage(e), "\n")
                NULL
            }
        )
        if (is.null(fit_model)) next

        # Extract batch ID from the filename to pass to the simulation function for naming and metadata purposes
        batch_id <- tryCatch(
            as.integer(sub(".*batch([0-9]+)_model\\.rds$", "\\1", basename(rds_path))),
            error = function(e) batch
        )

        # Run the simulation function
        simulation_results <- simulate_from_model(
            fit_obj  = fit_model,
            seed     = seed_k,
            batch_id = batch_id,
            family   = args$family
        )

        # Clean up the loaded model from memory to free up resources before the next iteration
        rm(fit_model); invisible(gc())

        # If simulation was successful, store the results in the list
        if (!is.null(simulation_results)) {
            replicat_results[[batch]] <- simulation_results
            cat(sprintf("    OK — %.1f min\n",
                        as.numeric(difftime(Sys.time(), t0, units = "min"))))
        } else {
            cat("    Batch", batch, "skipped.\n")
        }
    }

    # ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    #  Merge batches for this replicate
    # ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # Keep only the successful simulations and print error if no batch succeded
    valid_simulation <- Filter(function(r) !is.null(r), replicat_results)
    cat(sprintf("\n  Merging %d/%d successful batches...\n",
                length(valid_simulation), length(rds_files)))
    if (length(valid_simulation) == 0) {
        cat("  All batches failed — skipping replicate", rep_k, "\n")
        next
    }

    # Ensure that all batches have the same set of peaks (rows) before merging
    common_peaks  <- Reduce(intersect, lapply(valid_simulation, function(batch) rownames(batch$counts))) # retrieve the common peaks
    merged_counts <- do.call(cbind, lapply(valid_simulation, function(batch) batch$counts[common_peaks, , drop = FALSE]))
    merged_meta   <- do.call(rbind, lapply(valid_simulation, function(batch) batch$meta))

    # Prepare metadata for the replicate for h5ad export
    merged_meta$cell_type <- as.character(merged_meta$cell_type)  # convert the factor to character vector of all the cell types present
    merged_meta$replicate <- rep_k                                 # replicate index, useful when merging replicates downstream
    merged_meta$donor_id  <- as.character(merged_meta$donor_id)   # convert the factor to character vector containing all the donors used for the simulation (used for pseudobulk to normalise by the number of donors)
    merged_meta$seed      <- seed_k                                # seed used for this replicate (reproducibility)
    merged_meta$ADNC_cat  <- args$ADNC_cat                         # ADNC group of the cell, useful when merging different ADNC groups
    merged_meta$synthetic <- TRUE                                  # flag to distinguish synthetic cells from real ones

    cat(sprintf("  Merged: %d peaks x %d cells\n",
                nrow(merged_counts), ncol(merged_counts)))
    cat("  Donors in merged_meta :", paste(unique(merged_meta$donor_id), collapse = ", "), "\n")
    # ------------------------------------------------------------------------------------------------------------
    #  Save the merge values of the replicate as h5ad
    # ------------------------------------------------------------------------------------------------------------
    simulation_Anndata_file <- AnnData(
        X   = t(as.matrix(merged_counts)),
        obs = merged_meta,
        var = data.frame(peak = common_peaks, row.names = common_peaks)
    )

    #Create the name of the simulated file in format <ADNC>_rep<number>_<seed>
    out_name <- sprintf("%s_rep%02d_seed%d_scDesign3.h5ad",
                        args$ADNC_cat, rep_k + args$start_replicate_count - 1, seed_k)
    out_path <- file.path(args$outdir, out_name)
    simulation_Anndata_file$write_h5ad(out_path)
    cat("  Saved:", out_path, "\n")

    # Clean up merged results from memory before the next replicate to free up resources
    rm(simulation_Anndata_file, merged_counts, merged_meta, valid_simulation, replicat_results)
    invisible(gc())
}

# Print total elapsed time for the entire simulation process across all replicates
elapsed <- difftime(Sys.time(), start_time, units = "min")
cat(sprintf("\nDone — %d replicates in %.1f min\n", N_REP, as.numeric(elapsed)))
