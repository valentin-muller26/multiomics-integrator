library(ggplot2) meme chose12:30rcran_pkgs <- c(
  "purrr", "reshape2", "dplyr", "stringr", "psych", "ggplot2"
)

bioc_pkgs <- c(
  "DESeq2", "edgeR", "limma"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

new_cran <- cran_pkgs[!cran_pkgs %in% installed.packages()[, "Package"]]
if (length(new_cran)) {
  install.packages(new_cran)
}

new_bioc <- bioc_pkgs[!bioc_pkgs %in% installed.packages()[, "Package"]]
if (length(new_bioc)) {
  BiocManager::install(new_bioc, update = FALSE, ask = FALSE)
}
library(purrr)
library(DESeq2)
library(edgeR)
library(limma)
library(reshape2)
library(dplyr)
library(stringr)
library(psych)
library(ggplot2)



preprocess_rna <- function(rna_data, donors, star_qc_rows){
  #' Preprocess the raw RNA-seq count matrix
  #'
  #' @description
  #' Function that keeps only the columns corresponding to the donors of
  #' interest and removes the STAR QC rows (non-gene rows such as
  #' N_unmapped, N_multimapping, N_noFeature, N_ambiguous).
  #'
  #' @param rna_data wide format RNA-seq count dataframe. Rows = features
  #'   (ENSEMBL gene IDs + STAR QC rows), columns = samples HC<num><treatment>.
  #' @param donors character vector of donor IDs to keep (e.g. c("HC16","HC17","HC18","HC19")).
  #' @param star_qc_rows character vector of non-gene rownames to remove  (e.g. c("N_unmapped","N_multimapping","N_noFeature","N_ambiguous")).
  #' @return filtered RNA-seq dataframe with only donor columns and gene rows.
  list_donor <- paste(donors, collapse = "|")
  
  # Keep only donor columns
  rna_filtered <- rna_data %>%
    dplyr::select(matches(list_donor))
  
  # Remove STAR QC rows
  rna_filtered <- rna_filtered[!rownames(rna_filtered) %in% star_qc_rows, ]
  return(rna_filtered)
}



preprocess_atac <- function(atac_data, donors){
  #' Preprocess the raw ATAC-seq count matrix
  #'
  #' @description
  #' Function that cleans the column names so they match the RNA-seq sample naming convention HC<num><treatment>
  #' builds peak IDs of the form  chr<nb>_start_end as rownames, keeps only donor sample columns and
  #' removes peaks located on alternative contigs (rownames containing more  than two underscores)
  #' @param atac_data wide format ATAC-seq count dataframe
  #' @param donors character vector of donor IDs, ordered by replicate number (e.g. c("HC16","HC17","HC18","HC19") -> REP1 = HC16 ...) 
   #' @return atac_data wide format dataframe with column names in the format HC<num><treatment> and rownames = chr<nb>_start_end
  
  # Clean the name and generate the same treatment name as RNA-seq
  colnames(atac_data) <- gsub("\\.mLb\\.clN\\.sorted\\.bam", "", colnames(atac_data))
  colnames(atac_data) <- gsub("IF", "IFNg", colnames(atac_data))
  colnames(atac_data) <- gsub("IL413IFNg", "IL413IFN", colnames(atac_data))
  
  # Change the name to have the same format as the RNA-seq data
  # REP<i> is replaced by donors[i]
  for (i in seq_along(donors)) {
    pattern     <- paste0("^(.*)_REP", i, "$")
    replacement <- paste0(donors[i], "\\1")
    colnames(atac_data) <- gsub(pattern, replacement, colnames(atac_data))
  }
  
  # Use chromosome number _ start _ end as identifier
  atac_data <- atac_data %>%
    mutate(ID = paste0(Chr, "_", Start, "_", End))
  rownames(atac_data) <- atac_data$ID
  atac_data <- atac_data %>% dplyr::select(contains("HC"))
  
  # Remove alternative notation (peaks on alt contigs have extra underscores)
  n_underscores <- str_count(rownames(atac_data), "_")
  valid_peaks   <- rownames(atac_data)[n_underscores == 2]
  atac_data     <- atac_data[valid_peaks, , drop = FALSE]
  return(atac_data)
}


generate_long_dataframe <- function(dt, omic_name){
  #' Convert a wide omics matrix to MOFA2 long format
  #'
  #' @description
  #' Takes a wide matrix/dataframe (features in rows, samples in columns)
  #' and reshapes it to the long format expected by MOFA2.
  #' the feature name = to the gene ENSEMBL ID ou peak ID
  #' sample = sample name HC<number><treatment>
  #' view = the omics name (rna or atac)
  #' value = expression value of the RNA or ATAC seq feature
  #' group equal single group for all sample since did not use the multigroup option of MOFA due to having not enough samples
  #'
  #' @param dt Wide-format dataframe of the omic. Rows = features
  #'   (ENSEMBL gene IDs or peak IDs), columns = samples (HC<id><treatment>).
  #' @param omic_name Character. Name of the omic ("rna" or "atac").
  #' @return A long dataframe with columns: sample, group, feature, view, value.

  dt$feature <- rownames(dt)
  long_df <- melt(dt, id.vars = "feature",
                  variable.name = "sample",
                  value.name = "value")
  long_df$view <- omic_name
  # Single group: since did not use multi group option of MOFA 2 and donor effect already removed via limma
  long_df$group <- "single_group"
  # Reorder the column to have the order that MOFA attend 
  long_df <- long_df[, c("sample", "group", "feature", "view", "value")]
  return(long_df)
}

# ---------------------------------------------------------------------------
# Normalization only (no MAD yet)
# ---------------------------------------------------------------------------
normalise_rna <- function(dt){
  #' Normalise the RNA omics by applying Variance Stabilizing Transformation
  #'
  #' @description
  #' Function that filters genes with no expression in any sample and applies
  #' VST transformation on the dataframe.
  #'
  #' @param dt RNA-seq wide format dataframe. Rows = features
  #'   (ENSEMBL gene IDs), columns = samples (HC<id><treatment>).
  #' @return vst_dt VST transformed RNA-seq wide format matrix.
  dt <- dt[rowSums(dt) > 0, ]
  coldata <- data.frame(sample = colnames(dt), row.names = colnames(dt))
  dds_dt <- DESeqDataSetFromMatrix(dt, coldata, ~1)
  message("RNA: vst")
  vst_dt <- assay(vst(dds_dt, blind = TRUE))
  return(vst_dt)
}

normalise_atac <- function(dt, min_cpm = 1, min_samples = 3){
  #' Normalise the ATAC omics by TMM + logCPM transformation
  #'
  #' @description
  #' Function that filters peaks with no signal in any sample, removes
  #' low-count peaks (CPM >= min_cpm in at least min_samples samples),
  #' applies TMM normalisation and returns the logCPM matrix.
  #'
  #' @param dt ATAC-seq wide format dataframe. Rows = features
  #'   (peak IDs), columns = samples (HC<id><treatment>).
  #' @param min_cpm Minimum CPM threshold for the low-count filter.
  #' @param min_samples Minimum number of samples in which a peak must
  #'   reach min_cpm to be kept.
  #' @return logcpm TMM-normalised logCPM ATAC-seq wide format matrix.
  
  # Remove low expression peak 
  dt <- dt[rowSums(dt) > 0, ]
  dge <- DGEList(counts = dt)
  keep <- rowSums(cpm(dge) >= min_cpm) >= min_samples
  message(paste("ATAC: keeping", sum(keep), "/", length(keep), "peaks after low-count filter"))
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  
  # Normalise with TMM and logCPM
  message("ATAC: TMM")
  dge <- calcNormFactors(dge, method = "TMM")
  message("ATAC: logCPM")
  logcpm <- cpm(dge, log = TRUE, prior.count = 1)
  return(logcpm)
}

build_meta <- function(omic_data){
  #' Build a sample-level metadata dataframe from an omics matrix
  #'
  #' @description
  #' Function that receives a wide-format omics matrix/dataframe and extracts
  #' sample-level metadata from the column names:
  #' sample = donor id and treatment in format HC<num><treatment>
  #' donor = donor id extracted from the sample name (HC<num>)
  #' condition = treatment extracted from the sample name (<treatment>)
  #' rownames of the returned dataframe = sample names, so the metadata
  #' can be indexed directly by sample.
  #'
  #' @param omic_data wide format matrix or dataframe of the omics data.
  #'   Rows = features (ENSEMBL gene IDs or peak IDs in format
  #'   chr<nb>_start_end), columns = samples (HC<num><treatment>).
  #' @return dataframe containing the sample-level metadata of omic_data.
  sample_name = colnames(omic_data)
  data.frame(
    sample    = sample_name,
    donor     = str_extract(sample_name, "HC[0-9]+"),
    condition = str_remove(sample_name, "HC[0-9]+"),
    row.names = sample_name,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
#  Batch correction (limma) — preserves condition signal
# ---------------------------------------------------------------------------
remove_donor_effect <- function(omic_data, omic_name = ""){
   #' Remove donor effect from an omics matrix using limma
  #'
  #' @description
  #' Function that removes the donor (batch) effect from a normalised omics matrix 
  #' Builds the sample-level metadata from the column names, constructs
  #' a design matrix on the condition, and calls limma::removeBatchEffect with donor as the batch variable.
  #'
  #' @param omic_data wide format normalised omics matrix/dataframe.
  #'   Rows = features (ENSEMBL gene IDs or peak IDs), columns = samples (HC<num><treatment>).
  #' @param omic_name character, name of the omic ("rna" or "atac") used only for the log message.
  #' @return corrected matrix of the same dimensions as omic_data, with the donor effect removed.
  meta   <- build_meta(omic_data)
  design <- model.matrix(~ condition, data = meta)
  message(omic_name, ": removeBatchEffect on ", ncol(omic_data),
          " samples, ", length(unique(meta$donor)), " donors, ",
          length(unique(meta$condition)), " conditions")
  corrected <- removeBatchEffect(
    as.matrix(omic_data),
    batch  = meta$donor,
    design = design
  )
  return(corrected)
}


# ---------------------------------------------------------------------------
# MAD-based feature selection (AFTER batch correction)
# ---------------------------------------------------------------------------
select_top_mad <- function(omic_data, n_top = 10000, omic_name = ""){
  #' Select the top n features based on Median Absolute Deviation (MAD)
  #'
  #' @description
  #' Function that ranks features by their MAD across samples and keeps the n_top most variable ones. 
  #'If omic_data has fewer than n_top features, all features are returned.
  #'
  #' @param omic_data wide format normalised omics matrix/dataframe.
  #'   Rows = features (ENSEMBL gene IDs or peak IDs), columns = samples (HC<num><treatment>).
  #' @param n_top integer, number of features to keep.
  #' @param omic_name character, name of the omic ("rna" or "atac") used only for the log message.
  #' @return dataframe containing the n_top selected features.
  message(omic_name, ": selecting top ", n_top, " features by MAD")
  top_feature <- order(apply(omic_data, 1, mad), decreasing = TRUE)[1:min(n_top, nrow(omic_data))]
  return(as.data.frame(omic_data[top_feature, ]))
}

# ---------------------------------------------------------------------------
# PCA diagnostic
# ---------------------------------------------------------------------------
plot_pca_diagnostic <- function(omic_data, title = ""){
  #' PCA diagnostic plot of an omics matrix
  #'
  #' @description
  #' Function that performs a PCA on the samples of an omics matrix andplots PC1 vs PC2, 
  #' colors represent the donor and shape the condition
  #' the percentage of variance explaine by each PC is added to the axes
  #'
  #' @param omic_data wide format omics matrix/dataframe. Rows = features
  #'   (ENSEMBL gene IDs or peak IDs), columns = samples (HC<num><treatment>).
  #' @param title character, title of the plot.
  #' @return ggplot object of the PC1 vs PC2 scatter plot.
  meta <- build_meta(omic_data)
  pca  <- prcomp(t(as.matrix(omic_data)), center = TRUE, scale. = FALSE)
  var_pc <- round(100 * (pca$sdev^2) / sum(pca$sdev^2), 1)
  df <- data.frame(
    PC1       = pca$x[, 1],
    PC2       = pca$x[, 2],
    donor     = meta$donor,
    condition = meta$condition
  )
  
  # Manual shape
  shapes_8 <- c(16, 17, 15, 18, 8, 11, 10, 12)
  
  ggplot(df, aes(x = PC1, y = PC2, color = donor, shape = condition)) +
    geom_point(size = 3, alpha = 0.85) +
    scale_shape_manual(values = shapes_8) +
    labs(
      title = title,
      x = paste0("PC1 (", var_pc[1], "%)"),
      y = paste0("PC2 (", var_pc[2], "%)")
    ) +
    theme_bw()
}
