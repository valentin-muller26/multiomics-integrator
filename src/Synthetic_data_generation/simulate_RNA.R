library(Matrix)
library(dplyr)
library(zellkonverter)
library(SingleCellExperiment)
library(anndata)
library(argparse)
reticulate::py_install("anndata")
if (!require("scDesign3", quietly = TRUE))
    devtools::install_github("SONGDONGYUAN1994/scDesign3")
library(scDesign3)

# --------------------------------------------------------------------------- #
#  Parse arguments
# --------------------------------------------------------------------------- #
parser <- ArgumentParser(description = "scDesign3 RNA-seq — generate replicates from saved models")
parser$add_argument("-m", "--modeldir", type = "character", required = TRUE,
                    help = "Directory containing saved .rds model files")
parser$add_argument("-o", "--outdir",   type = "character", required = TRUE,
                    help = "Output directory for simulated h5ad replicates")
parser$add_argument("--donor_id",       type = "character", default = "donor",
                    help = "Donor/group identifier added to obs metadata [default: 'donor']")
parser$add_argument("--n_rep",          type = "integer",   default = 10,
                    help = "Number of replicates to generate [default: 10]")
args <- parser$parse_args()

N_REP      <- args$n_rep
SEED_START <- 100

dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)

start_time <- Sys.time()
cat("donor_id  :", args$donor_id, "\n")
cat("n_rep     :", N_REP, "\n\n")

# --------------------------------------------------------------------------- #
#  Discover model files
# --------------------------------------------------------------------------- #
manifest_path <- file.path(args$modeldir, "model_manifest.tsv")

if (file.exists(manifest_path)) {
    manifest  <- read.table(manifest_path, sep = "\t", header = TRUE,
                            stringsAsFactors = FALSE)
    rds_files <- manifest$rds_path
    cat("Loaded manifest:", manifest_path, "\n")
} else {
    rds_files <- sort(list.files(args$modeldir, pattern = "\\.rds$", full.names = TRUE))
    cat("No manifest found — using all .rds files in directory.\n")
}

cat("Found", length(rds_files), "model file(s).\n\n")
if (length(rds_files) == 0) stop("No model files found in: ", args$modeldir)

# --------------------------------------------------------------------------- #
#  Helper: simulate new counts from fitted marginals + copula
#
#  KEY DESIGN: We pass newCovariate exactly as construct_data() produced it.
#  simu_new() internally splits cells by corr_by groups (cell_type) and
#  relies on matching row structure between input_data and new_covariate.
#  Randomness between replicates comes from set.seed() before simu_new(),
#  which controls the Gaussian copula sampling.
# --------------------------------------------------------------------------- #
simulate_from_model <- function(fit_obj, seed, batch_id) {
    set.seed(seed)

    # Use the EXACT newCovariate from construct_data — no resampling
    new_cov <- fit_obj$dat$newCovariate

    cat("    new_cov dims       :", nrow(new_cov), "x", ncol(new_cov), "\n")
    if ("cell_type" %in% colnames(new_cov)) {
        ct_tab <- table(new_cov$cell_type)
        cat("    cell_type distrib  :",
            paste(names(ct_tab), ct_tab, sep = "=", collapse = ", "), "\n")
    }

    # Step 1 — extract parameters (deterministic given model + covariate)
    cat("    Extracting parameters...\n")
    para <- tryCatch(
        extract_para(
            sce           = fit_obj$sce_ref,
            marginal_list = fit_obj$marginal_fit,
            n_cores       = 1,
            family_use    = "nb",
            new_covariate = new_cov,
            data          = fit_obj$dat$dat
        ),
        error = function(e) {
            cat("    ERROR in extract_para:", conditionMessage(e), "\n")
            NULL
        }
    )
    if (is.null(para)) return(NULL)

    cat("    mean_mat dims      :", nrow(para$mean_mat), "x", ncol(para$mean_mat), "\n")

    # Step 2 — simulate counts (stochastic via copula, controlled by seed)
    cat("    Simulating counts...\n")
    counts <- tryCatch(
        simu_new(
            sce               = fit_obj$sce_ref,
            mean_mat          = para$mean_mat,
            sigma_mat         = para$sigma_mat,
            zero_mat          = para$zero_mat,
            quantile_mat      = NULL,
            copula_list       = fit_obj$copula_fit$copula_list,
            n_cores           = 1,
            family_use        = "nb",
            input_data        = fit_obj$dat$dat,
            new_covariate     = new_cov,
            important_feature = fit_obj$copula_fit$important_feature,
            filtered_gene     = fit_obj$dat$filtered_gene
        ),
        error = function(e) {
            cat("    ERROR in simu_new:", conditionMessage(e), "\n")
            NULL
        }
    )

    if (is.null(counts)) return(NULL)
    if (ncol(counts) == 0) {
        cat("    WARNING: simu_new returned 0 cells — skipping.\n")
        return(NULL)
    }

    cat("    simu_new output    :", nrow(counts), "genes x", ncol(counts), "cells\n")

    # Name cells uniquely
    actual_ncell     <- ncol(counts)
    ct_vec           <- as.character(new_cov$cell_type[seq_len(actual_ncell)])
    counter          <- ave(seq_along(ct_vec), ct_vec, FUN = seq_along)
    new_names        <- paste0(ct_vec, "_", counter, "_b", batch_id)
    colnames(counts) <- new_names

    # Build metadata
    meta           <- new_cov[seq_len(actual_ncell), , drop = FALSE]
    rownames(meta) <- new_names

    list(counts = counts, meta = meta)
}

# --------------------------------------------------------------------------- #
#  Generate replicates
# --------------------------------------------------------------------------- #
for (rep_k in seq_len(N_REP)) {
    seed_k <- SEED_START + rep_k
    cat(sprintf("\n======  Replicate %d/%d  (seed=%d)  ======\n",
                rep_k, N_REP, seed_k))

    rep_results <- vector("list", length(rds_files))

    for (j in seq_along(rds_files)) {
        rds_path <- rds_files[[j]]
        cat(sprintf("  [Batch %d/%d] Loading: %s\n", j, length(rds_files),
                    basename(rds_path)))
        t0 <- Sys.time()

        fit_j <- tryCatch(
            readRDS(rds_path),
            error = function(e) {
                cat("    ERROR loading model:", conditionMessage(e), "\n")
                NULL
            }
        )
        if (is.null(fit_j)) next

        batch_id <- tryCatch(
            as.integer(sub(".*batch([0-9]+)_model\\.rds$", "\\1", basename(rds_path))),
            error = function(e) j
        )

        sim_j <- simulate_from_model(
            fit_obj  = fit_j,
            seed     = seed_k,
            batch_id = batch_id
        )

        rm(fit_j); invisible(gc())

        if (!is.null(sim_j)) {
            rep_results[[j]] <- sim_j
            cat(sprintf("    OK — %.1f min\n",
                        as.numeric(difftime(Sys.time(), t0, units = "min"))))
        } else {
            cat("    Batch", j, "skipped.\n")
        }
    }

    # ------------------------------------------------------------------ #
    #  Merge batches for this replicate
    # ------------------------------------------------------------------ #
    valid <- Filter(Negate(is.null), rep_results)
    cat(sprintf("\n  Merging %d/%d successful batches...\n",
                length(valid), length(rds_files)))

    if (length(valid) == 0) {
        cat("  All batches failed — skipping replicate", rep_k, "\n")
        next
    }

    common_genes  <- Reduce(intersect, lapply(valid, function(r) rownames(r$counts)))
    merged_counts <- do.call(cbind, lapply(valid, function(r) r$counts[common_genes, , drop = FALSE]))
    merged_meta   <- do.call(rbind, lapply(valid, function(r) r$meta))

    # Clean up metadata for h5ad export
    merged_meta$cell_type <- as.character(merged_meta$cell_type)
    if ("donor_id" %in% colnames(merged_meta)) {
        merged_meta$donor_id <- as.character(merged_meta$donor_id)
    }
    merged_meta$replicate <- rep_k
    merged_meta$seed      <- seed_k
    merged_meta$donor_id  <- args$donor_id
    merged_meta$synthetic <- TRUE

    cat(sprintf("  Merged: %d genes x %d cells\n",
                nrow(merged_counts), ncol(merged_counts)))

    # ------------------------------------------------------------------ #
    #  Save as h5ad
    # ------------------------------------------------------------------ #
    sim_ad <- AnnData(
        X   = t(as.matrix(merged_counts)),
        obs = merged_meta,
        var = data.frame(gene = common_genes, row.names = common_genes)
    )

    out_name <- sprintf("%s_rep%02d_seed%d_scDesign3.h5ad",
                        args$donor_id, rep_k, seed_k)
    out_path <- file.path(args$outdir, out_name)
    sim_ad$write_h5ad(out_path)
    cat("  Saved:", out_path, "\n")

    rm(rep_results, valid, merged_counts, merged_meta, sim_ad)
    invisible(gc())
}

cat("\nAll replicates done in",
    round(difftime(Sys.time(), start_time, units = "min"), 2), "min\n")