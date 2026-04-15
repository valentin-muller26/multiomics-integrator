library(DESeq2)
library(ggplot2)
library(argparse)

parser <- ArgumentParser(description = "PCA analysis of pseudobulk data for QC of simulated vs real datasets")
parser$add_argument("-Sr", "--simRNA",  type = "character", required = TRUE,
                    help = "Path to simulated RNA pseudobulk matrix (tab-delimited, features x samples)")
parser$add_argument("-Sa", "--simATAC", type = "character", required = TRUE,
                    help = "Path to simulated ATAC pseudobulk matrix (tab-delimited, features x samples)")
parser$add_argument("-Rr", "--realRNA",  type = "character", required = TRUE,
                    help = "Path to real RNA pseudobulk matrix (tab-delimited, features x samples)")
parser$add_argument("-Ra", "--realATAC", type = "character", required = TRUE,
                    help = "Path to real ATAC pseudobulk matrix (tab-delimited, features x samples)")
parser$add_argument("-o", "--outdir", type = "character", required = TRUE,
                    help = "Output directory for PCA plots")
args <- parser$parse_args()

dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)

# Load data — check.names=FALSE to preserve _ in column names
real_rna  <- read.table(args$realRNA,  sep="\t", header=TRUE, row.names=1, check.names=FALSE)
sim_rna   <- read.table(args$simRNA,   sep="\t", header=TRUE, row.names=1, check.names=FALSE)
real_atac <- read.table(args$realATAC, sep="\t", header=TRUE, row.names=1, check.names=FALSE)
sim_atac  <- read.table(args$simATAC,  sep="\t", header=TRUE, row.names=1, check.names=FALSE)

vst_normalize <- function(real_df, sim_df) {
    # Find common features and filter out those with zero counts across all samples for the DESeq2 VST
    common <- intersect(rownames(real_df), rownames(sim_df))
    cat("Common features:", length(common), "\n")
    merged <- cbind(real_df[common, ], sim_df[common, ])
    merged <- round(as.matrix(merged))
    merged <- merged[rowSums(merged) > 0, ]
    cat("Features after zero-filtering:", nrow(merged), "\n")

    col_data <- data.frame(
        condition = factor(rep("A", ncol(merged))),
        row.names = colnames(merged)
    )
    dds <- DESeqDataSetFromMatrix(merged, colData=col_data, design=~1)
    assay(vst(dds, blind=TRUE))
}

parse_adnc <- function(sample_names) {
    # "real_Not_AD" -> "Not AD", "real_High" -> "High", "High_rep01" -> "High"
    adnc <- sub("^real_", "", sample_names)   # remove "real_" prefix
    adnc <- sub("_rep.*", "", adnc)           # remove "_rep..." suffix
    adnc <- gsub("Not_AD", "Not AD", adnc)    # fix Not_AD -> Not AD
    adnc
}

run_pca <- function(vst_mat, title, out_file) {
    # Perform PCA on the transposed VST matrix (samples as rows) and create a plot colored by ADNC levels
    pca <- prcomp(t(vst_mat), scale.=TRUE)
    var <- summary(pca)$importance[2, 1:2] * 100

    df_plot <- data.frame(
        PC1    = pca$x[, 1],
        PC2    = pca$x[, 2],
        sample = colnames(vst_mat),
        ADNC   = parse_adnc(colnames(vst_mat)),
        origin = ifelse(grepl("^real_", colnames(vst_mat)), "Real", "Simulated")
    )

    cat("ADNC levels found:", paste(unique(df_plot$ADNC), collapse=", "), "\n")

    p <- ggplot(df_plot, aes(x=PC1, y=PC2, color=ADNC)) +
        geom_point(data=subset(df_plot, origin=="Simulated"),
                   aes(shape=origin), size=2, alpha=0.5) +
        geom_point(data=subset(df_plot, origin=="Real"),
                   aes(shape=origin), size=4, stroke=0.8) +
        scale_color_manual(values=c(
            "Not AD"       = "#4C72B0",
            "Low"          = "#55A868",
            "Intermediate" = "#C44E52",
            "High"         = "#8172B2"
        )) +
        scale_shape_manual(values=c("Real"=16, "Simulated"=17)) +
        labs(title=title,
             x=sprintf("PC1 (%.1f%%)", var[1]),
             y=sprintf("PC2 (%.1f%%)", var[2])) +
        theme_bw()

    path_out <- file.path(args$outdir, out_file)
    ggsave(path_out, p, width=7, height=6, dpi=150)
    cat("Saved:", path_out, "\n")
}

cat("\n=== RNA ===\n")
vst_rna  <- vst_normalize(real_rna, sim_rna)
run_pca(vst_rna, "PCA — RNA VST (real vs simulated)", "PCA_RNA_vst.png")

cat("\n=== ATAC ===\n")
vst_atac <- vst_normalize(real_atac, sim_atac)
run_pca(vst_atac, "PCA — ATAC VST (real vs simulated)", "PCA_ATAC_vst.png")