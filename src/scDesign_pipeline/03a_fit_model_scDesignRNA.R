#' Script to fit scDesign3 models for RNA-seq data, batched by cell type
#'@param rnafile Path to the input h5ad RNA-seq file for the ADNC group
#'@param outdir Output directory to save the fitted models (.rds files)
#'@param batch_size Number of cell types to include in each batch (default: 5)
#'@param cells_per_type Number of cells to sample per cell type (default: 20)
#'@param seed Random seed for reproducibility (default: 42)
#'@param var_quantile Quantile for variance filtering (default: 0.95, i.e., keep top 95% most variable genes)
#'@examples     
#' Rscript 03a_fit_model_scDesignRNA.R \
#'    --rnafile "/data/users/vmuller/0_master_thesis/data/SEAD_Dataset/patient_subsets/group/Not AD.h5ad" \
#'    --outdir "/data/users/vmuller/0_master_thesis/data/SEAD_Dataset/models/RNA/Not AD"


library(Matrix)
library(dplyr)
library(zellkonverter)
library(SingleCellExperiment)
library(anndata)
library(argparse)
library(scDesign3)
 
# --------------------------------------------------------------------------- #
#  Parse arguments
# --------------------------------------------------------------------------- #
parser <- ArgumentParser(description = "scDesign3 RNA-seq — fit models (batched by cell type)")
parser$add_argument("-r", "--rnafile",  type = "character", required = TRUE,
                    help = "Path to the input h5ad RNA-seq file")
parser$add_argument("-o", "--outdir",   type = "character", required = TRUE,
                    help = "Output directory for saved models (.rds files)")
parser$add_argument("--batch_size", type = "integer", default = 5,
                    help = "Number of cell types to include in each batch (default: 5)")
parser$add_argument("--cells_per_type", type = "integer", default = 20,
                    help = "Number of cells to sample per cell type (default: 20)")
parser$add_argument("--seed", type = "integer", default = 42,
                    help = "Random seed for reproducibility (default: 42)")
parser$add_argument("--var_quantile", type = "double", default = 0.95,
                    help = "Quantile for variance filtering (default: 0.95, i.e., keep top 95% most variable genes)")

args <- parser$parse_args()
 

# --------------------------------------------------------------------------- 
#  Fit function
# --------------------------------------------------------------------------- 
fit_batch <- function(keep_types, rna_counts, cell_type, cell_type_list,
                      cells_per_type, seed) {
    #' Function to fit the scDesign3 model for a batch of cell types.
    #' This function subsets the data to the specified cell types, fit the model for the batch and return the model with the metadata about the batch.
    #' @param keep_types Vector of cell types to include in this batch (e.g., c("Microglia", "Astrocyte"))
    #' @param rna_counts Matrix of RNA counts (genes x cells)
    #' @param cell_type Named vector of cell types for each cell (names are cell barcodes)
    #' @param cell_type_list List of cell types, each containing the indices of the cells belonging to that type
    #' @param cells_per_type Number of cells to sample per cell type (e.g., 20)
    #' @param seed Random seed for reproducibility (e.g., 42)
    #' @return A list containing the fitted marginal and copula models, the reference SingleCellExperiment, the constructed data, and the barcodes of the cells used for fitting.
    #' @examples
    #' fit_result <- fit_batch(
    #'     keep_types = c("Microglia", "Astrocyte"),
    #'     rna_counts = rna_counts,
    #'     cell_type = cell_type,
    #'     cell_type_list = cell_type_list,
    #'     cells_per_type = 20,
    #'     seed = 42
    #' )    

    # Set seed for reproducibility
    set.seed(seed)
    
    # Subset the number of cells per type to a maximum of cells_per_type (e.g., 20) per donor to speed up the fitting process
    cells_keep <- unlist(lapply(keep_types, function(current_cell_type) { # external Loop over each cell type in the batch
        cell_indices <- cell_type_list[[current_cell_type]] # Retrieve the indices of cells belonging to the current cell type
        donors <- colData(sce_rna)$Donor.ID[cell_indices] # Get the donor IDs for those cells
 
        unlist(lapply(unique(donors), function(donor_name) { # internal loop over each donor for the current cell type
            donor_cell_indices <- cell_indices[donors == donor_name] # Retrieve the indices of cells belonging to the current donor
            n_take <- min(cells_per_type, length(donor_cell_indices)) # Determine how many cells to take from this donor 
            if (n_take == length(donor_cell_indices)) { # If the number of cells to take is equal to the number of available cells for this donor, take all of them
                donor_cell_indices
            } else {
                donor_cell_indices[sample.int(length(donor_cell_indices), n_take)] # Otherwise, randomly sample n_take cells from this donor
            }
        }))
    }))
    cells_keep    <- sort(unique(cells_keep))   # Sort and keep unique cell indices to avoid possible duplicates across all cell types in the batch
    rna_sub       <- rna_counts[, cells_keep, drop = FALSE]
    barcodes_keep <- colnames(rna_sub)

    cat("  Fitting on:", ncol(rna_sub), "cells\n")
 
    # metadata about the cells used for fitting (cell type and donor ID) -> needed for the model covariates
    cell_meta <- data.frame(
        cell_type = factor(cell_type[barcodes_keep]),
        donor_id  = factor(colData(sce_rna)$Donor.ID[cells_keep]),
        row.names = barcodes_keep
    )
    
    # Create SingleCellExperiment object for the subsetted data and the corresponding metadata for fitting the model
    rna_sce <- SingleCellExperiment(
        assays  = list(counts = as(rna_sub, "dgCMatrix")),
        colData = cell_meta
    )
 
    # Step 1 — construct data for model fitting (one row per cell, one column per gene, plus covariates)
    cat("  Constructing data...\n")
    dat <- construct_data(
        sce              = rna_sce,
        assay_use        = "counts",
        celltype         = "cell_type",
        pseudotime       = NULL,
        spatial          = NULL,
        other_covariates = "donor_id",
        corr_by          = "1"
    )
 
    # Step 2 — fit marginal distributions
    # Covariates for the mean (mu) formula: cell type + donor ID (to capture differences between donors)
    # Covariates for the dispersion (sigma) formula: 1 (constant dispersion across all cells) to reduce complexity and speed up fitting
    cat("  Fitting marginals...\n")
    marginal_fit <- fit_marginal(
        data            = dat,
        mu_formula      = "cell_type + donor_id",
        sigma_formula   = "1",
        family_use      = "nb",
        n_cores         = 1,
        usebam          = FALSE,
        parallelization = "pbmcmapply",
        trace           = FALSE
    )
 
    # Step 3 — fit Gaussian copula
    cat("  Fitting copula...\n")
    copula_fit <- fit_copula(
        sce           = rna_sce,
        assay_use     = "counts",
        marginal_list = marginal_fit,
        family_use    = "nb",
        copula        = "gaussian",
        n_cores       = 1,
        input_data    = dat$dat
    )

    # Return the fitted models (marginal_fit and copula_fit) and the metadata about the batch (the subsetted SingleCellExperiment, the constructed data, and the barcodes of the cells used for fitting) as a list
    list(    
        marginal_fit  = marginal_fit,
        copula_fit    = copula_fit,
        sce_ref       = rna_sce,
        dat           = dat,
        barcodes_keep = barcodes_keep
    )
}

# Define the main parameters for the fitting process
BATCH_SIZE     <- args$batch_size
CELLS_PER_TYPE <- args$cells_per_type
SEED           <- args$seed
VAR_QUANTILE   <- 1 - args$var_quantile

# Create output directory if it doesn't exist
dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)
 
# Time the entire fitting process
start_time <- Sys.time()
 
# --------------------------------------------------------------------------- 
#  Load data
# --------------------------------------------------------------------------- 
cat("Loading:", args$rnafile, "\n")
sce_rna    <- readH5AD(args$rnafile)
rna_counts <- assay(sce_rna, "UMIs") # Extract the RNA counts matrix (genes x cell)
cell_type  <- setNames(as.character(colData(sce_rna)$Subclass), colnames(sce_rna)) 
 
# Build per-cell-type index
cell_barcode   <- colnames(rna_counts) 
cell_type_list <- split(seq_along(cell_barcode), cell_type[cell_barcode])
# Cell type list ->A list of cell types, each containing the indices of the cells belonging to that type
# $Microglia
# [1]   3  17  42  88 ...
# $Astrocyte
# [1]   1   5  12  29 ...
# $Oligo
# [1]   2   4   9  15 ...
all_types      <- names(cell_type_list) # unique names of the cell types present in the dataset
cat("Found", length(all_types), "cell types:", paste(all_types, collapse = ", "), "\n\n")
 
# --------------------------------------------------------------------------- 
#   Variance filter (keep top 95% most variable genes) — computed on ALL cells
# --------------------------------------------------------------------------- 
set.seed(SEED)
gene_var    <- rowVars(rna_counts)   # sparseMatrixStats — no dense conversion
var_cutoff  <- quantile(gene_var, VAR_QUANTILE)
genes_kept  <- rownames(rna_counts)[gene_var >= var_cutoff]
rna_counts <- rna_counts[genes_kept, ]
cat("Genes after variance filter (top",
    paste0((1-VAR_QUANTILE) * 100, "% most variable):"),
    nrow(rna_counts), "\n\n")

    
# Save the kept gene names — needed in 04a_simulate_RNA.R to subset the simulated data to the same genes
genes_kept_path <- file.path(args$outdir, "genes_kept.rds")
saveRDS(genes_kept, genes_kept_path)
cat("Kept genes saved:", genes_kept_path, "\n\n")

# Subset sce_rna to the filtered genes so training_cells.h5ad is consistent
sce_rna <- sce_rna[genes_kept, ]

# Clean up gene_var to free memory before fitting the models
rm(gene_var); invisible(gc())

# Split into batches
batches   <- split(all_types, ceiling(seq_along(all_types) / BATCH_SIZE))
n_batches <- length(batches)
cat("Running", n_batches, "batches of up to", BATCH_SIZE, "cell types each.\n\n")

# --------------------------------------------------------------------------- #
#  Run all batches and save each model
# --------------------------------------------------------------------------- #
all_barcodes_used <- character(0)

for (i in seq_len(n_batches)) {
    celltypes_in_batch <- batches[[i]]
    cat(sprintf("\n=== Batch %d/%d  [%s] ===\n", i, n_batches, paste(celltypes_in_batch, collapse = ", ")))
    t0 <- Sys.time()
 
    fit_i <- tryCatch(
        fit_batch(
            keep_types     = celltypes_in_batch,
            rna_counts     = rna_counts,
            cell_type      = cell_type,
            cell_type_list = cell_type_list,
            cells_per_type = CELLS_PER_TYPE,
            seed           = SEED
        ),
        error = function(e) {
            cat("  ERROR:", conditionMessage(e), "\n")
            NULL
        }
    )
    
    if (!is.null(fit_i)) {
        # Add metadata about the batch to the fitted model list (cell types included, batch ID, number of cells per type, random seed, 
        #and the barcodes of the cells used for fitting)
        fit_i$.meta <- list(
            cell_types     = celltypes_in_batch,
            batch_id       = i,
            cells_per_type = CELLS_PER_TYPE,
            seed           = SEED,
            barcodes_keep  = fit_i$barcodes_keep
        )

        # Accumulate barcodes across all batches
        all_barcodes_used <- union(all_barcodes_used, fit_i$barcodes_keep)

        # Save the fitted model for this batch as an .rds file in the output directory with a name indicating the batch number (e.g., batch01_model.rds)
        rds_path <- file.path(args$outdir, sprintf("batch%02d_model.rds", i))
        saveRDS(fit_i, rds_path)
        cat(sprintf("  Model saved: %s\n", rds_path))
        cat(sprintf("  Done in %.1f min\n",
                    as.numeric(difftime(Sys.time(), t0, units = "min"))))
    } else {
        cat("  Batch", i, "skipped (error).\n")
    }
 
    rm(fit_i); invisible(gc())
}

# --------------------------------------------------------------------------- #
#  Export all training cells as a single h5ad
# --------------------------------------------------------------------------- #
cat("\nExporting training cells to h5ad...\n")
h5ad_path <- file.path(args$outdir, "training_cells.h5ad")
sce_train <- sce_rna[, all_barcodes_used]
writeH5AD(sce_train, h5ad_path)
cat("Training cells h5ad saved:", h5ad_path, "\n")
cat("Total training cells:", ncol(sce_train), "\n")

# --------------------------------------------------------------------------- #
#  Write manifest
# --------------------------------------------------------------------------- #
manifest <- data.frame(
    batch_id   = seq_len(n_batches),
    rds_path   = file.path(args$outdir, sprintf("batch%02d_model.rds", seq_len(n_batches))),
    cell_types = sapply(batches, paste, collapse = "|"),
    stringsAsFactors = FALSE
)
manifest_path <- file.path(args$outdir, "model_manifest.tsv")
write.table(manifest, manifest_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nManifest written:", manifest_path, "\n")
 
cat("\nTotal fitting time:",
    round(difftime(Sys.time(), start_time, units = "min"), 2), "min\n")
cat("Models saved in:", args$outdir, "\n")