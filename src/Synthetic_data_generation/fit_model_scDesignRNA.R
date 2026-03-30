library(Matrix)
library(dplyr)
library(zellkonverter)
library(SingleCellExperiment)
library(anndata)
library(argparse)
 
if (!require("scDesign3", quietly = TRUE))
    devtools::install_github("SONGDONGYUAN1994/scDesign3")
library(scDesign3)
 
# --------------------------------------------------------------------------- #
#  Parse arguments
# --------------------------------------------------------------------------- #
parser <- ArgumentParser(description = "scDesign3 RNA-seq — fit models (batched by cell type)")
parser$add_argument("-r", "--rnafile",  type = "character", required = TRUE,
                    help = "Path to the input h5ad RNA-seq file")
parser$add_argument("-o", "--outdir",   type = "character", required = TRUE,
                    help = "Output directory for saved models (.rds files)")
args <- parser$parse_args()
 
BATCH_SIZE     <- 5
CELLS_PER_TYPE <- 20
SEED           <- 42
 
dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)
 
start_time <- Sys.time()
 
# --------------------------------------------------------------------------- #
#  Load data
# --------------------------------------------------------------------------- #
sce_rna    <- readH5AD(args$rnafile)
rna_counts <- assay(sce_rna, "UMIs")
cell_type  <- setNames(as.character(colData(sce_rna)$Subclass), colnames(sce_rna))
 
# Build per-cell-type index
cell_barcode   <- colnames(rna_counts)
cell_type_list <- split(seq_along(cell_barcode), cell_type[cell_barcode])
all_types      <- names(cell_type_list)
cat("Found", length(all_types), "cell types:", paste(all_types, collapse = ", "), "\n\n")
 
# Split into batches
batches   <- split(all_types, ceiling(seq_along(all_types) / BATCH_SIZE))
n_batches <- length(batches)
cat("Running", n_batches, "batches of up to", BATCH_SIZE, "cell types each.\n\n")
 

fit_batch <- function(keep_types, rna_counts, cell_type, cell_type_list,
                      cells_per_type, seed) {
    set.seed(seed)
 
    cells_keep <- unlist(lapply(keep_types, function(current_cell_type) {
        cell_indices <- cell_type_list[[current_cell_type]]
        donors <- colData(sce_rna)$Donor.ID[cell_indices]
 
        unlist(lapply(unique(donors), function(donor_name) {
            donor_cell_indices <- cell_indices[donors == donor_name]
            n_take <- min(cells_per_type, length(donor_cell_indices))
            if (n_take == length(donor_cell_indices)) {
                donor_cell_indices
            } else {
                donor_cell_indices[sample.int(length(donor_cell_indices), n_take)]
            }
        }))
    }))
    cells_keep    <- sort(unique(cells_keep))
    rna_sub       <- rna_counts[, cells_keep, drop = FALSE]
    barcodes_keep <- colnames(rna_sub)

    cat("  Fitting on:", ncol(rna_sub), "cells\n")
 

    cell_meta <- data.frame(
        cell_type = factor(cell_type[barcodes_keep]),
        donor_id  = factor(colData(sce_rna)$Donor.ID[cells_keep]),
        row.names = barcodes_keep
    )
 
    rna_sce <- SingleCellExperiment(
        assays  = list(counts = as(rna_sub, "dgCMatrix")),
        colData = cell_meta
    )
 
    # Step 1 — construct covariate data
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
 
    # Step 2 — fit marginal distributions (one model per gene)
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
 
    list(
        marginal_fit  = marginal_fit,
        copula_fit    = copula_fit,
        sce_ref       = rna_sce,
        dat           = dat,
        barcodes_keep = barcodes_keep
    )
}
 
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
        fit_i$.meta <- list(
            cell_types     = celltypes_in_batch,
            batch_id       = i,
            cells_per_type = CELLS_PER_TYPE,
            seed           = SEED,
            barcodes_keep  = fit_i$barcodes_keep
        )

        # Accumulate barcodes across all batches
        all_barcodes_used <- union(all_barcodes_used, fit_i$barcodes_keep)

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