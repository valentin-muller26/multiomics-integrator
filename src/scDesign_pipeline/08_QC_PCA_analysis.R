#' PCA analysis of pseudobulk data for Quality Control of simulated vs real datasets
#'
#' Script that merge the real and simulated pseudobulk matrics and remove features with zero counts across all samples, 
#' then perform PCA and create a plot colored by ADNC levels and shape by origin (Real or Simulated).
#'
#' @param simRNA Path to simulated RNA pseudobulk matrix (tab-delimited, features x samples).
#' @param simATAC Path to simulated ATAC pseudobulk matrix (tab-delimited, features x samples).
#' @param realRNA Path to real RNA pseudobulk matrix (tab-delimited, features x samples).
#' @param realATAC Path to real ATAC pseudobulk matrix (tab-delimited, features x samples).
#' @param outdir Output directory for PCA plots.
#' @examples
#' Rscript 10_QC_PCA_analysis.R \
#'    --realRNA "$REAL_RNA" \
#'    --simRNA "$SIM_RNA" \
#'    --realATAC "$REAL_ATAC" \
#'    --simATAC "$SIM_ATAC" \
#'    --outdir "$OUTPUT"


suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(argparse)
})

## Load arguments from command line
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


merge_real_sim <- function(real_df, sim_df) {
    #' Function to merge real and simulated pseudobulk matrices, filter out features with zero counts, and return the merged matrix for PCA analysis.
    #' @param real_df Data frame of real pseudobulk matrix (features x samples).
    #' @param sim_df Data frame of simulated pseudobulk matrix (features x samples).
    #' @examples
    #' merged_matrix <- merge_real_sim(real_df, sim_df)

    # Find common features (genes or peaks) between real and simulated data
    common_features <- intersect(rownames(real_df), rownames(sim_df))
    cat("Common features:", length(common_features), "\n")
    cat("Features in real data:", nrow(real_df), "\n")
    cat("Features in simulated data:", nrow(sim_df), "\n")
    # Merge real and simulated data, round to integers
    merged <- cbind(real_df[common_features, ], sim_df[common_features, ])
    merged <- round(as.matrix(merged))
    # Filter out features with zero counts across all samples
    merged <- merged[rowSums(merged) > 0, ]
    cat("Features after zero-filtering:", nrow(merged), "\n")
    return(merged)
}

parse_adnc <- function(sample_names) {
    #' Function to parse ADNC levels from sample names and clean the names to return the correct format.
    #' @param sample_names Character vector of sample names (column names of the merged matrix)
    #' @examples
    #' adnc_levels <- parse_adnc(colnames(merged_matrix))
    #' example input: "real_Not_AD", "real_High", "High_rep01"
    #' desired output: "Not AD", "High", "High"
    adnc <- sub("^real_", "", sample_names)   # remove "real_" prefix
    adnc <- sub("_rep.*", "", adnc)           # remove "_rep..." suffix
    adnc <- gsub("Not_AD", "Not AD", adnc)    # fix Not_AD to  Not AD
    return(adnc)
}

run_pca <- function(merged_matrix, title, out_file) {
    #' Function to generate PCA plot from the merged pseudobulk matrix, colored by ADNC levels and shaped by origin (Real vs Simulated) 
    #' and save the plot to the output directory.
    #' @param merged_matrix Merged pseudobulk matrix (features x samples) after filtering.
    #' @param title Title for the PCA plot.
    #' @param out_file Filename for the output PCA plot (e.g., "PCA_RNA.png").
    #' @examples
    #' run_pca(merged_matrix, "PCA — RNA (real vs simulated)", "PCA_RNA.png")


    # Perform PCA on the transposed matrix (samples as rows) and create a plot colored by ADNC levels
    pca <- prcomp(t(merged_matrix), scale.=TRUE) 
    var <- summary(pca)$importance[2, 1:2] * 100

    df_plot <- data.frame(
        PC1    = pca$x[, 1],
        PC2    = pca$x[, 2],
        sample = colnames(merged_matrix),
        ADNC   = parse_adnc(colnames(merged_matrix)),
        origin = ifelse(grepl("^real_", colnames(merged_matrix)), "Real", "Simulated")
    )

    cat("ADNC levels found:", paste(unique(df_plot$ADNC), collapse=", "), "\n")

    p <- ggplot(df_plot, aes(x=PC1, y=PC2, color=ADNC)) +
        geom_point(data=subset(df_plot, origin=="Simulated"),
                   aes(shape=origin), size=2.5, alpha=0.5) +
        geom_point(data=subset(df_plot, origin=="Real"),
                   aes(shape=origin), size=5, stroke=0.8) +
        scale_color_manual(values=c(
            "Not AD"       = "#4C72B0",
            "Low"          = "#55A868",
            "Intermediate" = "#C44E52",
            "High"         = "#8172B2"
        )) +
        scale_shape_manual(values=c("Real"=16, "Simulated"=17)) +
        # Bigger legend keys/points so symbols stay readable when the figure is shrunk
        guides(
            color = guide_legend(override.aes = list(size = 5)),
            shape = guide_legend(override.aes = list(size = 5))
        ) +
        labs(title=title,
             x=sprintf("PC1 (%.1f%%)", var[1]),
             y=sprintf("PC2 (%.1f%%)", var[2])) +
        theme_bw(base_size = 18) +
        theme(
            plot.title   = element_text(size = 22, face = "bold"),
            axis.title   = element_text(size = 18),
            axis.text    = element_text(size = 18),
            legend.title = element_text(size = 20),
            legend.text  = element_text(size = 18),
            legend.key.size = unit(1.1, "lines")
        )

    #Save the PCA plot to the output directory
    path_out <- file.path(args$outdir, out_file)
    ggsave(path_out, p, width=8, height=6.5, dpi=300)
    cat("Saved:", path_out, "\n")
}


# Create output directory if it doesn't exist
dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)

# Load data — check.names=FALSE to preserve _ in column names
real_rna  <- read.table(args$realRNA,  sep="\t", header=TRUE, row.names=1, check.names=FALSE)
sim_rna   <- read.table(args$simRNA,   sep="\t", header=TRUE, row.names=1, check.names=FALSE)
real_atac <- read.table(args$realATAC, sep="\t", header=TRUE, row.names=1, check.names=FALSE)
sim_atac  <- read.table(args$simATAC,  sep="\t", header=TRUE, row.names=1, check.names=FALSE)

cat("\n=== RNA ===\n")
cat("Library sizes real:", "(min:", min(colSums(real_rna)), ", max:", max(colSums(real_rna)), ")\n")
cat("Library sizes sim:",  "(min:", min(colSums(sim_rna)), ", max:", max(colSums(sim_rna)), ")\n")
merged_rna  <- merge_real_sim(real_rna, sim_rna)
run_pca(merged_rna, "PCA : RNA (real vs simulated)", "PCA_RNA.png")

cat("\n=== ATAC ===\n")
cat("Library sizes real:", "(min:", min(colSums(real_atac)), ", max:", max(colSums(real_atac)), ")\n")
cat("Library sizes sim:",  "(min:", min(colSums(sim_atac)), ", max:", max(colSums(sim_atac)), ")\n")
merged_atac <- merge_real_sim(real_atac, sim_atac)
run_pca(merged_atac, "PCA : ATAC (real vs simulated)", "PCA_ATAC.png")