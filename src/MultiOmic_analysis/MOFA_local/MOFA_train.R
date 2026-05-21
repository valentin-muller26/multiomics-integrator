#---------------------------------------------------------------------------------------------------------------------------------------------------------
# Automatic installation of missing packages
#---------------------------------------------------------------------------------------------------------------------------------------------------------
# CRAN packages
cran_packages <- c("purrr", "reshape2", "dplyr", "stringr", "psych",
                   "ggplot2", "ggrepel", "openxlsx", "BiocManager")

# Bioconductor packages
bioc_packages <- c("MOFA2", "DESeq2", "edgeR", "limma",
                   "MOFAdata", "org.Hs.eg.db")

# Install missing CRAN packages
install_if_missing_cran <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing CRAN package: ", pkg)
    install.packages(pkg, dependencies = TRUE)
  }
}
invisible(sapply(cran_packages, install_if_missing_cran))

# Make sure BiocManager is available before installing Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install missing Bioconductor packages
install_if_missing_bioc <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing Bioconductor package: ", pkg)
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}
invisible(sapply(bioc_packages, install_if_missing_bioc))

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
if ("package:data.table" %in% search()) {
  detach("package:data.table", unload = TRUE, character.only = TRUE)
}

#MAIN VARIABLE for the filtering
star_qc_rows <- c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous")
donors <- c("HC16","HC17","HC18","HC19")


#Load the module containing the general function and function specific for mofa
source("utils.R")
source("MOFA_function.R")

# Set the main variable 
nb_factor <- 5
input_path_rna  <- "data/RNAseq_run3_AO1_HC14_HC19_geneCounts.txt"
input_path_atac <- "data/consensus_peaks.mLb.clN.featureCounts.txt"
outdir_model    <- "results/mofa/model"
outdir_graph    <- "results/mofa/graph"

#Create output directory if not present
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


# Run the model on 10 different seed and select the best on 
model_list =list()
seeds <- seq(1,10)
for (i in seq_along(seeds)) {
  model_list[[i]] <- train_mofa_model(omics_df, outdir_model, nb_factor = nb_factor, seed = seeds[i])
}
# Select best model and add metadata 
MOFAobject.trained <- select_model(model_list,plot=F)
MOFAobject.trained <- add_metadata(MOFAobject.trained)
print(compare_elbo(model_list))  

#Save the model and the dataframe for the RNA and ATAC seq
saveRDS(MOFAobject.trained, file.path(outdir_model, "MOFA_final_model.rds"))
saveRDS(rna,file.path(outdir_model, "rna.rds"))
saveRDS(atac,file.path(outdir_model, "atac.rds"))

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

print(elbo_plot)

ggsave(file.path(outdir_graph, "elbo_convergence.pdf"),
       plot = elbo_plot, width = 7, height = 5)


while (!is.null(dev.list())) dev.off()  # clean all the device 
# Comparison of the factor for the different model trained
pdf(file.path(outdir_graph, "compare_factors_stability.pdf"), width = 8, height = 7)
p <- compare_factors(model_list)
if (!is.null(p) && !is.null(p$gtable)) {
  grid::grid.draw(p$gtable)
} else {
  print(p)
}
while (!is.null(dev.list())) dev.off()  # clean all the device
