library(Matrix)
library(MatrixGenerics)
if (!require("sparseMatrixStats", quietly = TRUE)) {
    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("sparseMatrixStats", ask = FALSE, update = FALSE)
}
library(sparseMatrixStats)
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
parser <- ArgumentParser(description = "scDesign3 ATAC-seq — fit models (batched by cell type, multi-donor)")
parser$add_argument("-a", "--atacfile", type = "character", required = TRUE,
                    help = "Path to the input h5ad ATAC-seq file (multi-donor)")
parser$add_argument("-o", "--outdir",   type = "character", required = TRUE,
                    help = "Output directory for saved models (.rds files)")
args <- parser$parse_args()

BATCH_SIZE     <- 4
CELLS_PER_TYPE <- 20
SEED           <- 42
VAR_QUANTILE   <- 0.90

dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)

start_time <- Sys.time()

# --------------------------------------------------------------------------- #
#  Load data
# --------------------------------------------------------------------------- #
cat("Loading:", args$atacfile, "\n")
sce_atac    <- readH5AD(args$atacfile)
atac_counts <- assay(sce_atac, "X")
cell_type   <- setNames(as.character(colData(sce_atac)$Subclass), colnames(sce_atac))
donor_id    <- setNames(as.character(colData(sce_atac)$Donor.ID),  colnames(sce_atac))

cat("Loaded", nrow(atac_counts), "peaks x", ncol(atac_counts), "cells\n")
cat("Donors:", paste(unique(donor_id), collapse = ", "), "\n\n")

# Sort by cell type
col_names   <- colnames(atac_counts)
cell_order  <- order(cell_type[col_names])
atac_counts <- atac_counts[, cell_order]

# Build per-cell-type index
col_names      <- colnames(atac_counts)
cell_type_list <- split(seq_along(col_names), cell_type[col_names])
all_types      <- names(cell_type_list)
cat("Found", length(all_types), "cell types:", paste(all_types, collapse = ", "), "\n\n")

# --------------------------------------------------------------------------- #
#  Variance filter (top 10% most variable peaks) — computed on ALL cells
# --------------------------------------------------------------------------- #
set.seed(SEED)
peak_var    <- rowVars(atac_counts)   # sparseMatrixStats — no dense conversion
var_cutoff  <- quantile(peak_var, VAR_QUANTILE)
peaks_kept  <- rownames(atac_counts)[peak_var >= var_cutoff]
atac_counts <- atac_counts[peaks_kept, ]
cat("Peaks after variance filter (top", paste0((1 - VAR_QUANTILE) * 100, "%):"),
    nrow(atac_counts), "\n\n")

# Save the kept peak names — needed in simulate_replicates_atac.R
peaks_kept_path <- file.path(args$outdir, "peaks_kept.rds")
saveRDS(peaks_kept, peaks_kept_path)
cat("Kept peaks saved:", peaks_kept_path, "\n\n")

# Subset sce_atac to the filtered peaks so training_cells.h5ad is consistent
sce_atac <- sce_atac[peaks_kept, ]

rm(peak_var); invisible(gc())

# --------------------------------------------------------------------------- #
#  Split into batches
# --------------------------------------------------------------------------- #
batches   <- split(all_types, ceiling(seq_along(all_types) / BATCH_SIZE))
n_batches <- length(batches)
cat("Running", n_batches, "batches of up to", BATCH_SIZE, "cell types each.\n\n")

# --------------------------------------------------------------------------- #
#  Fit function
# --------------------------------------------------------------------------- #
fit_batch <- function(keep_types, atac_counts, cell_type, donor_id,
                      cell_type_list, cells_per_type, seed) {
    set.seed(seed)

    # Sample cells per (cell_type x donor)
    cells_keep <- unlist(lapply(keep_types, function(ct) {
        idx    <- cell_type_list[[ct]]
        donors <- donor_id[colnames(atac_counts)[idx]]

        unlist(lapply(unique(donors), function(d) {
            donor_idx <- idx[donors == d]
            n_take <- min(cells_per_type, length(donor_idx))
            if (n_take == length(donor_idx)) {
                donor_idx
            } else {
                donor_idx[sample.int(length(donor_idx), n_take)]
            }
        }))
    }))
    cells_keep <- sort(unique(cells_keep))

    atac_sub      <- atac_counts[, cells_keep, drop = FALSE]
    barcodes_keep <- colnames(atac_sub)
    cat("  Fitting on:", ncol(atac_sub), "cells\n")

    # cell_type and donor_id as factors (scDesign3 requirement)
    cell_meta <- data.frame(
        cell_type = factor(cell_type[barcodes_keep]),
        donor_id  = factor(donor_id[barcodes_keep]),
        row.names = barcodes_keep
    )

    atac_sce <- SingleCellExperiment(
        assays  = list(counts = as(atac_sub, "dgCMatrix")),
        colData = cell_meta
    )

    # Step 1 — construct covariate data
    cat("  Constructing data...\n")
    dat <- construct_data(
        sce              = atac_sce,
        assay_use        = "counts",
        celltype         = "cell_type",
        pseudotime       = NULL,
        spatial          = NULL,
        other_covariates = "donor_id",
        corr_by          = "1"
    )

    # Step 2 — fit marginal distributions (one model per peak)
    cat("  Fitting marginals...\n")
    marginal_fit <- fit_marginal(
        data            = dat,
        mu_formula      = "cell_type + donor_id",
        sigma_formula   = "1",
        family_use      = "nb",
        n_cores         = 1,
        parallelization = "pbmcmapply",
        trace           = TRUE
    )

    # Step 3 — fit copula
    cat("  Fitting copula\n")
    copula_fit <- fit_copula(
        sce           = atac_sce,
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
        sce_ref       = atac_sce,
        dat           = dat,
        barcodes_keep = barcodes_keep
    )
}

# --------------------------------------------------------------------------- #
#  Run all batches and save each model
# --------------------------------------------------------------------------- #
all_barcodes_used <- character(0)

for (i in seq_len(n_batches)) {
    bt <- batches[[i]]
    cat(sprintf("\n=== Batch %d/%d  [%s] ===\n", i, n_batches, paste(bt, collapse = ", ")))
    t0 <- Sys.time()

    fit_i <- tryCatch(
        fit_batch(
            keep_types     = bt,
            atac_counts    = atac_counts,
            cell_type      = cell_type,
            donor_id       = donor_id,
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
            cell_types     = bt,
            batch_id       = i,
            cells_per_type = CELLS_PER_TYPE,
            seed           = SEED
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
#  Export training cells as h5ad (filtered peaks + sampled cells only)
# --------------------------------------------------------------------------- #
cat("\nExporting training cells to h5ad...\n")
h5ad_path  <- file.path(args$outdir, "training_cells.h5ad")
sce_train  <- sce_atac[peaks_kept, all_barcodes_used]   # rows = filtered peaks
writeH5AD(sce_train, h5ad_path)
cat("Training cells h5ad saved:", h5ad_path, "\n")
cat("Peaks in training h5ad:  ", nrow(sce_train), "\n")
cat("Cells in training h5ad:  ", ncol(sce_train), "\n")

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