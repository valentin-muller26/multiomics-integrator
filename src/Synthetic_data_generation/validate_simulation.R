suppressPackageStartupMessages({
  library(argparse)
  library(SingleCellExperiment)
  library(scuttle)
  library(irlba)
  library(umap)
  library(CellMixS)
  library(MatrixGenerics)
  library(Rfast)
  library(reshape2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(cowplot)
  library(pbmcapply)
  library(zellkonverter)
  library(ggrastr)
  library(scales)
})

# =============================================================================
# 1. Arguments
# =============================================================================
parser <- ArgumentParser(description = "Validation of the simulated files")

parser$add_argument("--real_sce",
                    required = TRUE,
                    help     = "Path to the file used for fitting the model")
parser$add_argument("--sim_sce",
                    required = TRUE,
                    help     = "Path to the h5ad simulated filed")
parser$add_argument("--modality",
                    required = TRUE,
                    choices  = c("RNA", "ATAC"),
                    help     = "Modality RNA or ATAC")
parser$add_argument("--output_dir",
                    default  = ".",
                    help     = "Output directory")
parser$add_argument("--prefix",
                    default  = "validation",
                    help     = "Prefix of the output files")
parser$add_argument("--cell_type_col",
                    default  = "cell_type",
                    help     = "Name for the column data of the cell type")
parser$add_argument("--n_neighbors",
                    type     = "integer",
                    default  = 15,
                    help     = "n_neighbors for UMAP [default: 15]")
parser$add_argument("--min_dist",
                    type     = "double",
                    default  = 0.1,
                    help     = "min_dist for UMAP [default: 0.1]")
parser$add_argument("--n_pcs",
                    type     = "integer",
                    default  = 50,
                    help     = "Number of principal component [default: 50]")
parser$add_argument("--n_cores",
                    type     = "integer",
                    default  = 4,
                    help     = "Number of core use for the paralelisation [default: 4]")
parser$add_argument("--ngene_cor",
                    type     = "integer",
                    default  = 100,
                    help     = "Top N features pour la heatmap de correlation [defaut: 100]")
parser$add_argument("--normalization",
                    action   = "store_true",
                    help     = "Appliquer scuttle::logNormCounts avant PCA (RNA)")

args <- parser$parse_args()

# =============================================================================
# 2. Helpers
# =============================================================================

.load_sce <- function(path, label = path) {
  sce <- zellkonverter::readH5AD(path, use_hdf5 = FALSE)
  for (aname in c("UMIs", "counts", "X")) {
    if (aname %in% assayNames(sce)) {
      cat(sprintf("  Using assay '%s' for %s\n", aname, label))
      assayNames(sce)[assayNames(sce) == aname] <- "counts"
      break
    }
  }
  if (!"counts" %in% assayNames(sce))
    stop(sprintf("No count assay in %s. Available: %s",
                 label, paste(assayNames(sce), collapse = ", ")))
  if (!"cell_type" %in% colnames(colData(sce)) &&
       "Subclass"  %in% colnames(colData(sce))) {
    colData(sce)$cell_type <- colData(sce)$Subclass
    cat(sprintf("  Mapped 'Subclass' -> 'cell_type' in %s\n", label))
  }
  cat(sprintf("  %s : %d cells x %d features\n", label, ncol(sce), nrow(sce)))
  sce
}

# =============================================================================
# 3. Main function
# =============================================================================

validateSimulation <- function(real_sce,
                               sim_sce,
                               modality      = c("RNA", "ATAC"),
                               output_dir    = ".",
                               prefix        = "validation",
                               cell_type_col = "cell_type",
                               n_neighbors   = 15,
                               min_dist      = 0.1,
                               n_pcs         = 50,
                               n_cores       = 4,
                               ngene_cor     = 100,
                               normalization = FALSE) {

  modality      <- match.arg(modality)
  tfidf         <- (modality == "ATAC")
  sce_list_name <- c("real", "simulated")

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # -------------------------------------------------------------------------
  # 0. Cleaning and alignment
  # -------------------------------------------------------------------------
  message("=== 0. Cleaning ===")
  real_sce <- real_sce[, colSums(counts(real_sce)) > 0]
  sim_sce  <- sim_sce[,  colSums(counts(sim_sce))  > 0]

  common_features <- intersect(rownames(real_sce), rownames(sim_sce))
  if (length(common_features) == 0) stop("No common feature found.")
  message("  Common features : ", length(common_features))
  real_sce <- real_sce[common_features, ]
  sim_sce  <- sim_sce[common_features, ]

  sce_list <- list(real = real_sce, simulated = sim_sce)

  # -------------------------------------------------------------------------
  # 1. PCA
  # -------------------------------------------------------------------------
  message("=== 1. PCA ===")
  counts(real_sce) <- as.matrix(counts(real_sce))

  if (!tfidf) {
    if (normalization) {
      real_sce     <- scuttle::logNormCounts(real_sce, size.factors = NULL)
      test_pca_fit <- irlba::prcomp_irlba(t(logcounts(real_sce)),
                                           center = TRUE, scale. = FALSE, n = n_pcs)
    } else {
      test_pca_fit <- irlba::prcomp_irlba(t(log1p(counts(real_sce))),
                                           center = TRUE, scale. = FALSE, n = n_pcs)
    }
  } else {
    test_pca_fit <- irlba::prcomp_irlba(
      t(log1p(Signac::RunTFIDF(counts(real_sce)))),
      center = TRUE, scale. = TRUE, n = n_pcs)
  }

  sce_list <- lapply(sce_list, function(x) {
    counts(x) <- as.matrix(counts(x))
    if (!tfidf) {
      if (normalization) {
        x <- scuttle::logNormCounts(x, size.factors = NULL)
      } else {
        logcounts(x) <- log1p(counts(x))
      }
      reducedDim(x, "PCA") <- predict(test_pca_fit, newdata = t(logcounts(x)))
    } else {
      logcounts(x) <- log1p(counts(x))
      reducedDim(x, "PCA") <- predict(test_pca_fit,
                                       newdata = t(log1p(Signac::RunTFIDF(counts(x)))))
    }
    x
  })
  message("  PCA OK (", ncol(test_pca_fit$x), " PCs)")

  # -------------------------------------------------------------------------
  # 2. UMAP
  # -------------------------------------------------------------------------
  message("=== 2. UMAP ===")
  test_umap_fit <- umap::umap(test_pca_fit$x,
                               n_neighbors = n_neighbors,
                               min_dist    = min_dist)

  sce_list <- pbmcapply::pbmclapply(sce_list, function(x) {
    res <- predict(object = test_umap_fit, data = reducedDim(x, "PCA"))
    colnames(res) <- c("UMAP1", "UMAP2")
    reducedDim(x, "UMAP") <- res
    x
  }, mc.cores = n_cores)

  sce_list <- lapply(seq_along(sce_list), function(i) {
    sce <- sce_list[[i]]
    colData(sce)$Method <- sce_list_name[i]
    sce
  })
  names(sce_list) <- sce_list_name
  message("  UMAP OK")

  rd_list <- lapply(sce_list, function(x) {
    rd      <- colData(x) %>% as_tibble()
    rd_pca  <- reducedDim(x, "PCA")  %>% as_tibble()
    rd_umap <- reducedDim(x, "UMAP") %>% as_tibble()
    rd %>% bind_cols(rd_pca) %>% bind_cols(rd_umap)
  })
  rd_tbl <- bind_rows(rd_list, .id = "Method")
  message("  rd_tbl OK")

  # -------------------------------------------------------------------------
  # 4. CellMixS
  # -------------------------------------------------------------------------
  message("=== 4. CellMixS ===")
  ref_sce <- sce_list[["real"]]
  colData(ref_sce)$Method <- "real2"

  metric_res <- pbmcapply::pbmclapply(sce_list, function(x) {
    colnames(x) <- paste0("Cell", seq_len(ncol(x)))
    count_combine <- cbind(counts(x), counts(ref_sce))
    sce_combine   <- SingleCellExperiment(list(counts = count_combine))
    reducedDim(sce_combine, "PCA")  <- rbind(reducedDim(x, "PCA"),
                                              reducedDim(ref_sce, "PCA"))
    reducedDim(sce_combine, "UMAP") <- rbind(reducedDim(x, "UMAP"),
                                              reducedDim(ref_sce, "UMAP"))
    colData(sce_combine)$Method <- c(colData(x)$Method, colData(ref_sce)$Method)

    sce_pca  <- CellMixS::evalIntegration(
      metrics  = c("isi", "entropy"), sce_combine,
      k = 50, n_dim = 2, cell_min = 4,
      res_name = c("weighted_isi", "entropy"),
      group = "Method", dim_red = "PCA")
    sce_umap <- CellMixS::evalIntegration(
      metrics  = c("isi", "entropy"), sce_combine,
      k = 50, n_dim = 2, cell_min = 4,
      res_name = c("weighted_isi", "entropy"),
      group = "Method", dim_red = "UMAP")

    c(PCA_lisi     = mean(colData(sce_pca)$weighted_isi,  na.rm = TRUE),
      PCA_entropy  = mean(colData(sce_pca)$entropy,       na.rm = TRUE),
      UMAP_lisi    = mean(colData(sce_umap)$weighted_isi, na.rm = TRUE),
      UMAP_entropy = mean(colData(sce_umap)$entropy,      na.rm = TRUE))
  }, mc.cores = n_cores)

  metric_res <- simplify2array(metric_res)
  metric_tbl <- t(metric_res) %>%
    as_tibble() %>%
    mutate(Method = sce_list_name)

  message("  CellMixS OK")
  print(metric_tbl)

  # -------------------------------------------------------------------------
  # 5. Eight statistic metrics
  # -------------------------------------------------------------------------
  message("=== 5. Eight statistic metrics ===")

  mean_log_expr <- pbmcapply::pbmclapply(sce_list, function(x)
    MatrixGenerics::rowMeans2(as.matrix(logcounts(x))),
    mc.cores = n_cores)
  message("  mean_log_expr OK")

  var_log_expr <- pbmcapply::pbmclapply(sce_list, function(x)
    MatrixGenerics::rowVars(as.matrix(logcounts(x))),
    mc.cores = n_cores)
  message("  var_log_expr OK")

  gene_detect_freq <- pbmcapply::pbmclapply(sce_list, function(x)
    MatrixGenerics::rowMeans2(as.matrix(counts(x)) == 0),
    mc.cores = n_cores)
  message("  gene_detect_freq OK")

  gene_cor <- pbmcapply::pbmclapply(sce_list, function(x) {
    cor_mat <- Rfast::cora(as.matrix(logcounts(x)))
    as.vector(cor_mat[upper.tri(cor_mat)])
  }, mc.cores = n_cores)
  message("  gene_cor OK")

  log_lib_size <- pbmcapply::pbmclapply(sce_list, function(x)
    log(Rfast::colsums(as.matrix(counts(x)))),
    mc.cores = n_cores)
  message("  log_lib_size OK")

  cell_distance <- pbmcapply::pbmclapply(sce_list, function(x)
    c(Rfast::Dist(reducedDim(x, "PCA"), square = FALSE)),
    mc.cores = n_cores)
  message("  cell_distance OK")

  cell_detect_freq <- pbmcapply::pbmclapply(sce_list, function(x)
    MatrixGenerics::colMeans2(as.matrix(counts(x)) == 0),
    mc.cores = n_cores)
  message("  cell_detect_freq OK")

  cell_cor <- pbmcapply::pbmclapply(sce_list, function(x) {
    cor_mat <- Rfast::cora(t(as.matrix(logcounts(x))))
    as.vector(cor_mat[upper.tri(cor_mat)])
  }, mc.cores = n_cores)
  message("  cell_cor OK")

  metric_cols <- c("Mean log expression", "Var log expression",
                   "Gene detection freq", "Gene correlation",
                   "Log library size",    "Cell distance",
                   "Cell detection freq", "Cell correlation")

  metric_tbl <- metric_tbl %>%
    mutate(
      `Mean log expression` = mean_log_expr,
      `Var log expression`  = var_log_expr,
      `Gene detection freq` = gene_detect_freq,
      `Gene correlation`    = gene_cor,
      `Log library size`    = log_lib_size,
      `Cell distance`       = cell_distance,
      `Cell detection freq` = cell_detect_freq,
      `Cell correlation`    = cell_cor
    )
  message("  Eight statistic metrics OK")

  # -------------------------------------------------------------------------
  # 6. Heatmap
  # -------------------------------------------------------------------------
  message("=== 6. Heatmap ===")
  train_logcounts <- t(as.matrix(logcounts(sce_list[["real"]])))
  ngene_cor       <- min(ngene_cor, ncol(train_logcounts))

  # Filtrer les features a variance nulle avant de selectionner le top N
  feat_vars   <- apply(train_logcounts, 2, var)
  nonzero_var <- which(feat_vars > 0)
  if (length(nonzero_var) == 0) stop("All features have zero variance.")

  top_n         <- min(ngene_cor, length(nonzero_var))
  feature_order <- order(colSums(train_logcounts[, nonzero_var]), decreasing = TRUE)
  top_idx       <- nonzero_var[feature_order[seq_len(top_n)]]
  top_mask      <- seq_len(ncol(train_logcounts)) %in% top_idx
  message("  Top features with non-zero variance : ", sum(top_mask))

  corr_mat_ref <- cor(train_logcounts[, top_mask], use = "pairwise.complete")
  corr_mat_ref[is.na(corr_mat_ref)] <- 0   # securite : NA residuels -> 0
  hclust_res   <- hclust(dist(corr_mat_ref))

  corr_mat_list <- pbmcapply::pbmclapply(sce_list, function(x) {
    mat     <- t(as.matrix(logcounts(x)[top_mask, ]))
    cor_mat <- Rfast::cora(mat)
    cor_mat[is.na(cor_mat)] <- 0
    cor_mat[hclust_res$order, hclust_res$order]
  }, mc.cores = n_cores)
  names(corr_mat_list) <- sce_list_name

  cor_melted <- lapply(corr_mat_list, reshape2::melt)
  corr_dat   <- Reduce(rbind, cor_melted)
  corr_dat$Method <- Reduce(c, lapply(sce_list_name, function(nm)
    rep(nm, nrow(cor_melted[[nm]]))))
  corr_dat$Method <- factor(corr_dat$Method, levels = sce_list_name)
  message("  Heatmap OK")

  # -------------------------------------------------------------------------
  # 7. Figures
  # -------------------------------------------------------------------------
  message("=== 7. Figures ===")

  method_levels <- c("real", "simulated")
  method_labels <- c("Real data", "Simulated")
  pal <- setNames(scales::hue_pal()(2), method_labels)

  .relabel <- function(x)
    factor(x, levels = method_levels, labels = method_labels)

  .ks_d <- function(a, b) {
    if (length(a) == 0 || length(b) == 0) return(NA_real_)
    suppressWarnings(ks.test(a, b)$statistic)
  }

  ks_vals <- c(
    `Mean log expression` = .ks_d(mean_log_expr[["real"]],    mean_log_expr[["simulated"]]),
    `Var log expression`  = .ks_d(var_log_expr[["real"]],     var_log_expr[["simulated"]]),
    `Gene detection freq` = .ks_d(gene_detect_freq[["real"]], gene_detect_freq[["simulated"]]),
    `Gene correlation`    = .ks_d(gene_cor[["real"]],         gene_cor[["simulated"]]),
    `Log library size`    = .ks_d(log_lib_size[["real"]],     log_lib_size[["simulated"]]),
    `Cell distance`       = .ks_d(cell_distance[["real"]],    cell_distance[["simulated"]]),
    `Cell detection freq` = .ks_d(cell_detect_freq[["real"]], cell_detect_freq[["simulated"]]),
    `Cell correlation`    = .ks_d(cell_cor[["real"]],         cell_cor[["simulated"]])
  )

  message("  KS values :")
  print(ks_vals)

  # --- tbl_summarystats : construction robuste ---
  # (evite le crash vec_rbind quand unnest rencontre un vecteur vide)
  tbl_summarystats <- do.call(rbind, lapply(metric_cols, function(mc) {
    do.call(rbind, lapply(sce_list_name, function(meth) {
      idx  <- match(meth, metric_tbl$Method)
      vals <- metric_tbl[[mc]][[idx]]
      if (is.null(vals) || length(vals) == 0) return(NULL)
      tibble(
        Method = .relabel(meth),
        Metric = factor(mc, levels = metric_cols),
        Value  = as.numeric(vals)
      )
    }))
  }))
  tbl_summarystats <- tbl_summarystats %>%
    filter(!is.na(Metric), !is.na(Value)) %>%
    droplevels()

  message("  Niveaux Metric : ", nlevels(tbl_summarystats$Metric))

  # --- ks_tbl_long : une ligne par métrique, côté Simulated uniquement ---
  ypos_tbl <- tbl_summarystats %>%
    group_by(Metric) %>%
    summarise(
      ypos = if (all(is.na(Value))) NA_real_
             else max(Value, na.rm = TRUE) + 0.18 * diff(range(Value, na.rm = TRUE)),
      .groups = "drop"
    )

  ks_tbl_long <- tibble(
    Metric = factor(names(ks_vals), levels = metric_cols),
    ks_raw = unlist(ks_vals)
  ) %>%
    mutate(ks_label = case_when(
      is.na(ks_raw)  ~ NA_character_,
      ks_raw < 0.005 ~ "<0.01",
      TRUE           ~ format(round(ks_raw, 2), nsmall = 2)
    )) %>%
    left_join(ypos_tbl, by = "Metric") %>%
    filter(!is.na(Metric), !is.na(ks_label)) %>%
    droplevels() %>%
    mutate(Method = factor("Simulated", levels = method_labels))

  # --- p_summarystats ---
  p_summarystats <- ggplot(tbl_summarystats,
                            aes(x = Method, y = Value, color = Method)) +
    geom_violin(scale = "width", trim = TRUE, lwd = 1.1) +
    facet_wrap(~Metric, scales = "free_y", nrow = 2, ncol = 4) +
    theme_bw() +
    theme(aspect.ratio     = 1.5,
          legend.position  = "none",
          axis.text.x      = element_text(angle = 45, hjust = 1),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()) +
    scale_color_manual(values = pal) +
    xlab("") + ylab("") +
    geom_text(data        = ks_tbl_long,
              aes(x = Method, y = ypos, label = ks_label),
              angle = 45, vjust = 0, hjust = 0.5, size = 3,
              inherit.aes = FALSE, na.rm = TRUE)

  # --- corr_p ---
  dat_corr_list <- corr_dat %>%
    mutate(Method = .relabel(as.character(Method))) %>%
    group_split(Method)

  r_corr <- sapply(dat_corr_list, function(x) {
    if (nrow(dat_corr_list[[1]]) > 0 && nrow(x) > 0)
      cor(dat_corr_list[[1]]$value, x$value, method = "pearson")
    else
      NA_real_
  })

  dat_text_corr <- tibble(
    r      = r_corr,
    Method = .relabel(method_levels)
  ) %>%
    mutate(label = if_else(!is.na(r), paste0("r=", round(r, digits = 2)), "")) %>%
    mutate(label = if_else(as.character(Method) == method_labels[1], "", label))

  corr_p <- corr_dat %>%
    mutate(Method = .relabel(as.character(Method))) %>%
    ggplot(aes(x = Var2, y = Var1, fill = value)) +
    facet_wrap(~Method, nrow = 1) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                          midpoint = 0, limit = c(-1, 1), space = "Lab",
                          name = "Pearson\nCorrelation") +
    theme_bw() +
    theme(legend.position  = "bottom",
          axis.text        = element_blank(),
          axis.ticks       = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()) +
    xlab("") + ylab("") + coord_fixed() +
    geom_text(data        = dat_text_corr %>% filter(label != ""),
              aes(x = 100, y = 0, label = label),
              vjust = -0.5, hjust = 1.05, inherit.aes = FALSE)

  # --- p_PCA ---
  dat_text_PCA <- metric_tbl %>%
    mutate(label  = paste0("mLISI=", round(PCA_lisi, digits = 2)),
           Method = .relabel(Method)) %>%
    mutate(label  = if_else(as.character(Method) == method_labels[1], "", label)) %>%
    dplyr::select(Method, label) %>%
    filter(label != "")

  p_PCA <- rd_tbl %>%
    mutate(Method      = .relabel(Method),
           `Cell type` = .data[[cell_type_col]]) %>%
    ggplot(aes(x = PC1, y = PC2)) +
    ggrastr::rasterize(
      geom_point(size = 0.5, alpha = 0.5, aes(color = `Cell type`)), dpi = 300) +
    facet_wrap(~Method, nrow = 1) +
    theme_bw() +
    theme(legend.position  = "bottom",
          legend.title     = element_text(angle = 0),
          aspect.ratio     = 1,
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text        = element_blank(),
          axis.ticks       = element_blank()) +
    geom_text(data        = dat_text_PCA,
              aes(x = Inf, y = -Inf, label = label),
              vjust = -0.5, hjust = 1.05, inherit.aes = FALSE) +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))

  # --- p_UMAP ---
  dat_text_UMAP <- metric_tbl %>%
    mutate(label  = paste0("mLISI=", round(UMAP_lisi, digits = 2)),
           Method = .relabel(Method)) %>%
    mutate(label  = if_else(as.character(Method) == method_labels[1], "", label)) %>%
    dplyr::select(Method, label) %>%
    filter(label != "")

  p_UMAP <- rd_tbl %>%
    mutate(Method      = .relabel(Method),
           `Cell type` = .data[[cell_type_col]]) %>%
    ggplot(aes(x = UMAP1, y = UMAP2)) +
    ggrastr::rasterize(
      geom_point(size = 0.5, alpha = 0.5, aes(color = `Cell type`)), dpi = 300) +
    facet_wrap(~Method, nrow = 1) +
    theme_bw() +
    theme(legend.position  = "bottom",
          legend.title     = element_text(angle = 0),
          aspect.ratio     = 1,
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text        = element_blank(),
          axis.ticks       = element_blank()) +
    geom_text(data        = dat_text_UMAP,
              aes(x = Inf, y = -Inf, label = label),
              vjust = -0.5, hjust = 1.05, inherit.aes = FALSE) +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))

  p_final <- cowplot::plot_grid(
    p_summarystats, corr_p, p_PCA, p_UMAP,
    labels      = c("a", "b", "c", "d"),
    label_size  = 15,
    ncol        = 2,
    rel_heights = c(1, 1),
    align       = "hv",
    axis        = "lr"
  )

  # -------------------------------------------------------------------------
  # 8. Saving
  # -------------------------------------------------------------------------
  message("=== 8. Saving files ===")
  base_path <- file.path(output_dir, paste0(prefix, "_", modality))

  cowplot::ggsave2(paste0(base_path, "_panel.pdf"),
                   plot = p_final, width = 13, height = 10)
  ggsave(paste0(base_path, "_summarystats.png"),
         plot = p_summarystats, width = 12, height = 6,  dpi = 300)
  ggsave(paste0(base_path, "_corr_heatmap.png"),
         plot = corr_p,         width = 8,  height = 5,  dpi = 300)
  ggsave(paste0(base_path, "_pca.png"),
         plot = p_PCA,          width = 8,  height = 5,  dpi = 300)
  ggsave(paste0(base_path, "_umap.png"),
         plot = p_UMAP,         width = 8,  height = 5,  dpi = 300)

  readr::write_tsv(
    metric_tbl %>% dplyr::select(Method, PCA_lisi, PCA_entropy, UMAP_lisi, UMAP_entropy),
    paste0(base_path, "_cellmixs.tsv"))
  readr::write_tsv(rd_tbl,   paste0(base_path, "_rd.tsv"))
  readr::write_tsv(corr_dat, paste0(base_path, "_corr_dat.tsv"))

  message("  Files saved in : ", output_dir)

  invisible(list(
    rd_tbl     = rd_tbl,
    metric_tbl = metric_tbl,
    corr_dat   = corr_dat,
    sce_list   = sce_list,
    plots = list(
      summarystats = p_summarystats,
      corr         = corr_p,
      pca          = p_PCA,
      umap         = p_UMAP,
      panel        = p_final
    )
  ))
}

# =============================================================================
# Run
# =============================================================================
message("Loading the data")
message("  Real    : ", args$real_sce)
message("  Simule  : ", args$sim_sce)
message("  Modalite: ", args$modality)

real_sce <- .load_sce(args$real_sce, "real")
sim_sce  <- .load_sce(args$sim_sce,  "simulated")

validateSimulation(
  real_sce      = real_sce,
  sim_sce       = sim_sce,
  modality      = args$modality,
  output_dir    = args$output_dir,
  prefix        = args$prefix,
  cell_type_col = args$cell_type_col,
  n_neighbors   = args$n_neighbors,
  min_dist      = args$min_dist,
  n_pcs         = args$n_pcs,
  n_cores       = args$n_cores,
  ngene_cor     = args$ngene_cor,
  normalization = args$normalization
)

message("=== Finished ===")