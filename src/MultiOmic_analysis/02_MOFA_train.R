#' MOFA2 Training Pipeline for Multi-Omics Integration (RNA-seq + ATAC-seq)
#'
#' @description
#' Trains a MOFA2 model on paired RNA-seq and ATAC-seq data
#' The pipelinenormalizes each omic, removes the donor effect, selects the top variable
#' features using the mean absolute deviation (MAD) to have 10,000 feature in each omics , 
#' and fits MOFA2 across 10 random seeds to select the most stable model based on the ELBO criterion.
#'@param nb_factor Maximum number of MOFA2 latent factors to fit (default: 5)
#'@param input_path_rna Path to the raw RNA-seq count matrix (features x samples)
#'@param input_path_atac Path to the raw ATAC-seq consensus peak count matrix (features x samples)
#'@param outdir_model Output directory for the trained model and preprocessed matrices (.rds)
#'@param outdir_graph Output directory for the diagnostic plots (PDF)
#'@output MOFA_final_model.rds Best-ELBO MOFA2 model with sample metadata attached
#'@output rna.rds Final preprocessed RNA-seq matrix fed to MOFA2
#'@output atac.rds Final preprocessed ATAC-seq matrix fed to MOFA2
#'@output data_overview.pdf Overview of the MOFA input data matrix (samples × views)
#'@output qc_boxplots_pca.pdf Boxplots of the top-MAD RNA and ATAC matrices
#'@output pca_rna_before.pdf PCA of RNA before donor correction
#'@output pca_rna_after.pdf PCA of RNA after donor correction
#'@output pca_atac_before.pdf PCA of ATAC before donor correction
#'@output pca_atac_after.pdf PCA of ATAC after donor correction
#'@output elbo_convergence.pdf ELBO trajectory of the best MOFA2 model
#'@output compare_elbo_models.pdf ELBO comparison across the 10 trained models
#'@output compare_factors_stability.pdf Factor correlation heatmap across seeds
#'@examples
#' Rscript train_mofa.R \
#'    --input_path_rna "/data/users/vmuller/0_master_thesis/data/Davos/rna_counts.tsv" \
#'    --input_path_atac "/data/users/vmuller/0_master_thesis/data/Davos/atac_counts.tsv" \
#'    --outdir_model "/data/users/vmuller/0_master_thesis/results/MOFA/model" \
#'    --outdir_graph "/data/users/vmuller/0_master_thesis/results/MOFA/graphs" \
#'    --nb_factor 5

library(MOFA2)
library(purrr)
library(DESeq2)
library(edgeR)
library(limma)
library(reshape2)
library(dplyr)
library(stringr)
library(psych)
library(ggplot2)
library(MOFAdata)
library(ggrepel)
library(openxlsx)
library(org.Hs.eg.db)
library(argparse)
if ("package:data.table" %in% search()) {
  detach("package:data.table", unload = TRUE, character.only = TRUE)
}

parser <- ArgumentParser(description = "Train a MOFA2 model on multi-omics data (RNA-seq + ATAC-seq)")
parser$add_argument("-f", "--nb_factor",
                    type = "integer", default = 5,
                    help = "Maximum number of latent factors to fit. (default: 5)")

parser$add_argument("-om", "--outdir_model",
                    type = "character", required = TRUE,
                    help = "Output directory where the trained MOFA model (.hdf5) will be saved.")

parser$add_argument("-og", "--outdir_graph",
                    type = "character", required = TRUE,
                    help = "Output directory for plots generated during the training (elbo and comparison models).")

parser$add_argument("-ir", "--input_path_rna",
                    type = "character", required = TRUE,
                    help = "Path to the raw RNA-seq count matrix (features x samples).")

parser$add_argument("-ia", "--input_path_atac",
                    type = "character", required = TRUE,
                    help = "Path to the raw ATAC-seq consensus peak count matrix (features x samples)")

args <- parser$parse_args()

# Source the functions for MOSim simulation and data processing
path_scripts <- Sys.getenv("PATH_SCRIPTS")
if (path_scripts == "") {
    stop("Environment variable PATH_SCRIPTS is not set. Source 00_config.sh before running.")
}
source(file.path(path_scripts, "utils.R"))
source(file.path(path_scripts, "MOFA_function.R"))


#MAIN VARIABLE for the filtering
star_qc_rows <- c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous")
donors <- c("HC16","HC17","HC18","HC19")


# Set the main variable 
nb_factor <- args$nb_factor
outdir_model    <- args$outdir_model
outdir_graph <- args$outdir_graph
input_path_rna = args$input_path_rna
input_path_atac= args$input_path_atac

#Create the output directory if not present
dir.create(outdir_model, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_graph, recursive = TRUE, showWarnings = FALSE)
#---------------------------------------------------------------------------------------------------------------------------------------------------------
# RNA-seq pipeline: load -> preprocess -> normalize -> batch correct -> MAD
#---------------------------------------------------------------------------------------------------------------------------------------------------------
rna_raw <- read.table(input_path_rna,
                      header = TRUE, row.names = 1)
rna_raw <- preprocess_rna(rna_raw,donors,star_qc_rows)

# 1. Normalize (VST on all genes)
rna_norm <- normalise_rna(rna_raw)

# 2. Remove donor effect (on all genes, not just the top MAD ones)
rna_corrected <- remove_donor_effect(rna_norm, "RNA")

# 3. Select top variable features AFTER correction
rna <- select_top_mad(rna_corrected, n_top = 10000, omic_name = "RNA")


#---------------------------------------------------------------------------------------------------------------------------------------------------------
# ATAC-seq pipeline: load -> preprocess -> normalize -> batch correct -> MAD
#---------------------------------------------------------------------------------------------------------------------------------------------------------
atac_data <- read.table(input_path_atac,
                        header = TRUE, row.names = 1)
atac_data <- preprocess_atac(atac_data,donors)

# 1. Normalize (TMM + log-CPM on all peaks passing the cpm filter)
atac_norm <- normalise_atac(atac_data)

# 2. Remove donor effect
atac_corrected <- remove_donor_effect(atac_norm, "ATAC")

# 3. Select top variable features AFTER correction
atac <- select_top_mad(atac_corrected, n_top = 10000, omic_name = "ATAC")


#----------------------------------------------------------------------------------------------------------------------------------------------------------
# Sanity check + PCA diagnostic before/after
#----------------------------------------------------------------------------------------------------------------------------------------------------------
message("Common samples: ", length(intersect(colnames(rna), colnames(atac))))

# Boxplots on the final top-10k matrices
pdf(file.path(outdir_graph, "qc_boxplots_pca.pdf"))
par(mfrow = c(2, 1))
boxplot(rna,  las = 2, main = "RNA after batch correction + top MAD",  outline = FALSE)
boxplot(atac, las = 2, main = "ATAC after batch correction + top MAD", outline = FALSE)
par(mfrow = c(1, 1))
dev.off()

# PCA diagnostic on the same top-10k subset, before vs after correction
# We re-apply MAD on rna_norm/atac_norm just for a fair before/after comparison
rna_top_before  <- select_top_mad(rna_norm,  n_top = 10000, omic_name = "RNA (pre-corr)")
atac_top_before <- select_top_mad(atac_norm, n_top = 10000, omic_name = "ATAC (pre-corr)")

# Generate the different PCA
pca_rna_before  <- plot_pca_diagnostic(rna_top_before,  "RNA — before donor correction (top MAD)")
pca_rna_after   <- plot_pca_diagnostic(rna,             "RNA — after donor correction (top MAD)")
pca_atac_before <- plot_pca_diagnostic(atac_top_before, "ATAC — before donor correction (top MAD)")
pca_atac_after  <- plot_pca_diagnostic(atac,            "ATAC — after donor correction (top MAD)")

# Save the plots
ggsave(file.path(outdir_graph, "pca_rna_before.pdf"),  pca_rna_before,  width = 7, height = 5)
ggsave(file.path(outdir_graph, "pca_rna_after.pdf"),   pca_rna_after,   width = 7, height = 5)
ggsave(file.path(outdir_graph, "pca_atac_before.pdf"), pca_atac_before, width = 7, height = 5)
ggsave(file.path(outdir_graph, "pca_atac_after.pdf"),  pca_atac_after,  width = 7, height = 5)

#----------------------------------------------------------------------------------------------------------------------------------------------------------
# Long format and MOFA
#----------------------------------------------------------------------------------------------------------------------------------------------------------
long_rna  <- generate_long_dataframe(rna,  "rna")
long_atac <- generate_long_dataframe(atac, "atac")
omics_df  <- rbind(long_rna, long_atac)

# Save the data overview plot
MOFAobject_overview <- create_mofa(omics_df)
ggsave(
  file.path(outdir_graph, "data_overview.pdf"),
  plot_data_overview(MOFAobject_overview),
  width = 7, height = 4
)
rm(MOFAobject_overview)

# Run the model on 10 different seed and select the best on 
model_list =list()
seeds <- seq(1,10)
for (i in seq_along(seeds)) {
  model_list[[i]] <- train_mofa_model(omics_df, outdir_model, nb_factor = nb_factor, seed = seeds[i])
}
# Select best model and add metadata 
MOFAobject.trained <- select_model(model_list,plot=F)
MOFAobject.trained <- add_metadata(MOFAobject.trained)

#Save the model and the dataframe for the RNA and ATAC seq
saveRDS(MOFAobject.trained, file.path(outdir_model, "MOFA_final_model.rds"))
saveRDS(rna,file.path(outdir_model, "rna.rds"))
saveRDS(atac,file.path(outdir_model, "atac.rds"))

# Save ELBO comparison across the 10 trained models
pdf(file.path(outdir_graph, "compare_elbo_models.pdf"), width = 7, height = 5)
print(compare_elbo(model_list))
dev.off()

# Create the elbow for the list of model train and print the one for the best model
elbos <- sapply(model_list, get_elbo)
print(elbos)

# Difference best and other elbo 
elbos_centered <- elbos - max(elbos)
print(elbos_centered)

#Plot the elbo
elbo_df <- data.frame(
  iteration = seq_along(MOFAobject.trained@training_stats$elbo),
  elbo = MOFAobject.trained@training_stats$elbo
) |> filter(!is.nan(elbo))

elbo_plot <- ggplot(elbo_df, aes(x = iteration, y = elbo)) +
  geom_line(color = "steelblue") +
  labs(title = "MOFA2 ELBO convergence",
       x = "Iteration", y = "ELBO") +
  theme_bw()



ggsave(file.path(outdir_graph, "elbo_convergence.pdf"),
       plot = elbo_plot, width = 7, height = 5)


# Comparison of the factor stability across the 10 trained models
while (!is.null(dev.list())) dev.off()
pdf(file.path(outdir_graph, "compare_factors_stability.pdf"), width = 8, height = 7)
p <- compare_factors(model_list)
if (!is.null(p) && !is.null(p$gtable)) {
  grid::grid.draw(p$gtable)
} else {
  print(p)
}
while (!is.null(dev.list())) dev.off()
