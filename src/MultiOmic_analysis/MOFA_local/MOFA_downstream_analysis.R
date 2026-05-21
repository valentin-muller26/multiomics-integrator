
# Downstream analysis of a trained MOFA2 model (RNA-seq + ATAC-seq)
# This script performs post-training analysis of a MOFA2 model
#   1. Load the trained model and pre-computed RNA / ATAC data
#   2. Visualise explained variance and correlations with covariates
#   3. Explore latent factors (1D and 2D) by condition / donor
#   4. Plot heatmaps of top genes / peaks for condition-associated factors
#   5. Run gene set enrichment analysis (GSEA) on the RNA view
#   6. Inspect specific genes of interest in the factor space
#   7. Annotate ATAC-seq peaks with genomic features via ChIPseeker
 

# -------------------------------------------------------------------------------------
#  Package installation and loading
# -------------------------------------------------------------------------------------

cran_pkgs <- c(
  "purrr", "reshape2", "dplyr", "stringr", "psych",
  "ggplot2", "ggrepel", "openxlsx", "pheatmap", "tidyr"
)

bioc_pkgs <- c(
  "MOFA2", "MOFAdata", "DESeq2", "edgeR", "limma",
  "org.Hs.eg.db", "ChIPseeker", "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "GenomicRanges"
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

library(MOFA2)
library(MOFAdata)
library(DESeq2)
library(edgeR)
library(limma)
library(purrr)
library(dplyr)
library(tidyr)
library(reshape2)
library(stringr)
library(psych)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(openxlsx)
library(org.Hs.eg.db)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)

if ("package:data.table" %in% search()) {
  detach("package:data.table", unload = TRUE, character.only = TRUE)
}
source("utils.R")
source("MOFA_function.R")

# Set the output directory
outdir_model    <- "results/mofa/model"
outdir_graph    <- "results/mofa/graph"
output_file  <- file.path(outdir_graph, "pathway_enrich.xlsx")

#Create output directory if not present
dir.create(outdir_model, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_graph, recursive = TRUE, showWarnings = FALSE)

# Load the model and the dataset for RNA and ATAC
MOFAobject.trained <- readRDS(file.path(outdir_model, "MOFA_final_model.rds"))
rna  <- readRDS(file.path(outdir_model, "rna.rds"))
atac <- readRDS(file.path(outdir_model, "atac.rds"))


#General parameter for the plots
PCA_comparison    <- c("Factor1", "Factor3")
genes_of_interest <- list("ENSG00000184557", "ENSG00000185338") #list of the gene of interest for the Specific gene anaylsis


#----------------------------------------------------------------------------------------------------------------------------------------------------------
# MOFA plots
#----------------------------------------------------------------------------------------------------------------------------------------------------------
# Number of latent factors retained by the trained model
nb_factor_train   <- MOFAobject.trained@dimensions$K
# Names of all available factors (e.g. "Factor1", "Factor2", ...)
available_factors <- factors_names(MOFAobject.trained)


#---------------------------------------------------------------------------------------------------------------------------------------------
#Plots variance explained and correlation variance-covariate(treatment /donor)
#-------------------------------------------------------------------------------------------------------------------------------------------
 
# Variance explained by each factor in each omics view (RNA, ATAC)
p_var_view  <- plot_variance_explained(MOFAobject.trained, x = "view", y = "factor")
# Total variance explained per view (cumulative across factors)
p_var_total <- plot_variance_explained(MOFAobject.trained, plot_total = TRUE)[[2]]


# Extract the variance explained 
var_explained <- get_variance_explained(MOFAobject.trained)

# Extract the variance explained per factor and view(omics) (matrix :factor = line and view = column)
var_per_factor <- var_explained$r2_per_factor[[1]]
print(var_per_factor)

# Variance total 
var_total <- var_explained$r2_total[[1]]
print(var_total)

ggsave(file.path(outdir_graph, "variance_explained_per_factor.pdf"), p_var_view,  width = 6, height = 5)
ggsave(file.path(outdir_graph, "variance_explained_total.pdf"),     p_var_total, width = 6, height = 4)

#plot correlation between factor
pdf(file.path(outdir_graph, "factor_correlation.pdf"), width = 6, height = 5)
p_factor_cor <- plot_factor_cor(MOFAobject.trained)
if (inherits(p_factor_cor, "pheatmap")) {
  grid::grid.draw(p_factor_cor$gtable)
} else {
  print(p_factor_cor)
}
dev.off()

#plot correlation factor and covariate (donor, treatment)
p_cov <- correlate_factors_with_covariates(
  MOFAobject.trained,
  covariates = c("condition", "donor"),
  plot = "log_pval"
)
pdf(file.path(outdir_graph, "factor_covariate_association.pdf"), width = 5, height = 5)
print(p_cov)
dev.off()

cor_results <- correlate_factors_with_covariates(
  MOFAobject.trained,
  covariates = c("condition", "donor"),
  plot = "log_pval", return_data = TRUE
)


# Retrive the list of factor that are significantly linkewd to the treatement for plotting the heatmap for those factor
threshold <- -log10(0.05)
sig_factors_names     <- rownames(cor_results)[cor_results[, "condition"] > threshold]
list_factor_condition <- as.numeric(gsub("Factor", "", sig_factors_names))
sig_values            <- cor_results[sig_factors_names, "condition"]
print(list_factor_condition)
message("Factors significantly associated with condition (p.adj < 0.05):\n",
        paste0("  ", sig_factors_names, ": -log10(p.adj) = ", round(sig_values, 3),
               collapse = "\n"))






#----------------------------------------------------------------------------------------------------------------------------------------------------
# Factor visualisation (1D per factor, 2D scatter)
#----------------------------------------------------------------------------------------------------------------------------------------------------

# 1D plot: each factor on the x-axis, samples coloured by condition, shaped by donor.
# used to see which factor seperates which condition
p_factor_condition <- plot_factor(MOFAobject.trained,
                                  factors    = 1:nb_factor_train,
                                  color_by   = "condition",
                                  shape_by = "donor",
                                  dot_size   = 3,
                                  dodge      = TRUE,
                                  add_violin = FALSE)


# Get the factor values and the metadata
Z <- get_factors(MOFAobject.trained, factors = "all")[[1]]
meta <- samples_metadata(MOFAobject.trained)

# Merge factor values with metadata, then reshape to long format  
df <- as.data.frame(Z) %>%
  mutate(sample = rownames(.)) %>%
  left_join(meta, by = "sample") %>%
  pivot_longer(cols = starts_with("Factor"),
               names_to = "factor", values_to = "value")

# Show the value of the factor per condition
df %>%
  filter(factor %in% sig_factors_names) %>%
  group_by(factor, condition) %>%
  summarise(mean_value = mean(value), .groups = "drop") %>%
  pivot_wider(names_from = condition, values_from = mean_value)

#2D scatter plot of the two factors

p_factor_2d <- plot_factors(MOFAobject.trained,
                            factors  = PCA_comparison,
                            color_by = "condition",
                            shape_by = "donor",
                            dot_size = 3)

ggsave(file.path(outdir_graph, "factors_by_condition.pdf"), p_factor_condition, width = 10, height = 5)
ggsave(file.path(outdir_graph, paste0("factors_2D_", paste(PCA_comparison, collapse = "_vs_"), ".pdf")),
       p_factor_2d, width = 6, height = 5)




#--------------------------------------------------------------------------------------------------------------------------------------------
#Per-factor linked to treatement heatmaps and top feature weights
#---------------------------------------------------------------------------------------------------------------------------------------------
#initilised list of top gene for latter diagnosis ans ATAC peak annotation step
top_peaks <- list()
top_gene <- list()
for (factor in list_factor_condition) {
  factor_name <- paste0("Factor", factor)
  

  #Weight plot of the most contributing gene/peak for the current factor
  pdf(file.path(outdir_graph, paste0("factor", factor, "_weights_heatmaps.pdf")),
      width = 8, height = 6)
  
  print(plot_top_weights(MOFAobject.trained, view = "rna",  factor = factor_name, nfeatures = 20))
  print(plot_top_weights(MOFAobject.trained, view = "atac", factor = factor_name, nfeatures = 20))
  
  # heatmap
  hm_rna <- plot_data_heatmap(MOFAobject.trained,
                    view               = "rna",
                    factor             = factor_name,
                    features           = 25,
                    cluster_rows       = TRUE,
                    cluster_cols       = T,
                    annotation_samples = c("condition", "donor"),
                    main               = paste("Heatmap RNA:", factor_name))
  top_gene[[factor_name]] <- hm_rna$tree_row$labels[hm_rna$tree_row$order]
  hm_atac <-plot_data_heatmap(MOFAobject.trained,
                    view               = "atac",
                    factor             = factor_name,
                    features           = 25,
                    cluster_rows       = TRUE,
                    cluster_cols       = T,
                    annotation_samples = c("condition", "donor"),
                    main               = paste("Heatmap ATAC:", factor_name))
  #Retrieve the peak in the heatmap
  top_peaks[[factor_name]] <- hm_atac$tree_row$labels[hm_atac$tree_row$order]
  
  
  print(plot_data_scatter(MOFAobject.trained,
                          view       = "rna",
                          factor     = factor_name,
                          features   = 6,         # 6 top features RNA
                          color_by   = "condition",
                          shape_by   = "donor",
                          add_lm     = TRUE,    
                          dot_size   = 2.5))
  
  print(plot_data_scatter(MOFAobject.trained,
                          view       = "atac",
                          factor     = factor_name,
                          features   = 6,         # 6 top peaks ATAC
                          color_by   = "condition",
                          shape_by   = "donor",
                          add_lm     = TRUE,
                          dot_size   = 2.5))
  dev.off()
}

# Diagnostic of the correlation of the top gene between the factor
intersect(top_gene$Factor1,top_gene$Factor3)
intersect(top_gene$Factor1,top_gene$Factor5)
intersect(top_gene$Factor3,top_gene$Factor5)
W <- get_weights(MOFAobject.trained, views = "rna")[[1]]
cor(W[, c("Factor1", "Factor3", "Factor5")])



#---------------------------------------------------------------------------------------------------------------------------------------------------------
#  Gene Set Enrichment Analysis (GSEA)
#---------------------------------------------------------------------------------------------------------------------------------------------------------

#load the feature set bundled with MOFAdata
#   - reactomeGS : Reactome pathway database
#   - MSigDB C2  : curated gene sets
#   - MSigDB C5  : Gene Ontology gene sets
data("reactomeGS")
data("MSigDB_v6.0_C2_human")
data("MSigDB_v6.0_C5_human")
#Drop pathways with empty or missing names (cause errors in plotting)
reactomeGS_clean <- reactomeGS[!is.na(rownames(reactomeGS)) &
                                 rownames(reactomeGS) != "", ]

featuresets <- list(
  C2       = MSigDB_v6.0_C2_human,
  C5       = MSigDB_v6.0_C5_human,
  Reactome = reactomeGS_clean
)
signs <- c("positive", "negative", "all")

#vector containing containing the list of analysis resulting of no significant pathway (for diagnosis)
no_significant_pathways <- character()


excel_file <- createWorkbook()

for (set_name in names(featuresets)) {
  set <- featuresets[[set_name]]
  
  pdf(file.path(outdir_graph, paste0("GSEA_", set_name, ".pdf")),
      width = 14, height = 7)
  
  for (sign in signs) {
    enrichment <- run_enrichment(
      MOFAobject.trained,
      view             = "rna",
      factors          = list_factor_condition,
      feature.sets     = set,
      sign             = sign,
      statistical.test = "parametric"
    )
    tryCatch({
      p_hm <- plot_enrichment_heatmap(enrichment,
                                      alpha        = 0.05,
                                      max.pathways = 25) +
        ggtitle(paste0(set_name, " - ", sign, " weights"))
      print(p_hm)
    }, error = function(e) {
      message("No heatmap for : ", set_name, " - ", sign, " (", e$message, ")")
    })
    
    for (factor in sig_factors_names) {
      title <- paste0(factor, " - ", sign, " weights - ", set_name)
      
      tryCatch({
        p1 <- plot_enrichment(enrichment, factor = factor, max.pathways = 15) + ggtitle(title)
        print(p1)
        p2 <- plot_enrichment_detailed(enrichment, factor = factor, max.genes = 8, max.pathways = 15) + ggtitle(title)
        print(p2)
      }, error = function(e) {
        message("No significant pathway for : ", title)
        no_significant_pathways <<- c(no_significant_pathways, paste(title, "[plot_enrichment]"))
      })
      
      tryCatch({
        pathway_table <- extract_pathway_genes(enrichment, factor)
        
        if (!is.null(pathway_table) && nrow(pathway_table) > 0) {
          sheet_name <- substr(paste(set_name, sign, factor, sep = "_"), 1, 31)
          addWorksheet(excel_file, sheet_name)
          writeData(excel_file, sheet_name, pathway_table)
          setColWidths(excel_file, sheet_name, cols = 1:ncol(pathway_table),
                       widths = c(50, 12, 12, 14, 14, 80))
          freezePane(excel_file, sheet_name, firstRow = TRUE)
        }
      }, error = function(e) {
        message("Failed to extract pathway genes for : ", title, " — ", e$message)
      })
    }
  }
  
  dev.off()
}

print(no_significant_pathways)
saveWorkbook(excel_file, output_file, overwrite = TRUE)
message("Excel file save : ", output_file)



#----------------------------------------------------------------------------------------------------------------------
#Specific gene anaylsis
#-------------------------------------------------------------------------------------------------------------------

p_gene_expr <- genes_of_interest %>% map(~ plot_factors(MOFAobject.trained,
                                            factors  = PCA_comparison,
                                            color_by = .,
                                            scale    = TRUE,
                                            legend   = TRUE) +
                               geom_text_repel(aes(label = sample),
                                               size = 3,
                                               max.overlaps = Inf,
                                               min.segment.length = 0)
) %>% cowplot::plot_grid(plotlist = ., nrow = 1)

ggsave(file.path(outdir_graph, "factors_2D_by_gene_expression.pdf"),
       p_gene_expr, width = 12, height = 5)



#----------------------------------------------------------------------------------------------------------------------------------
#Atac seq anotation
#----------------------------------------------------------------------------------------------------------------------------------
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

wb <- createWorkbook()

for (factor_name in names(top_peaks)) {
  top_peak_factor <- top_peaks[[factor_name]]
  
  # Conversion to Grange for chipseeker
  parts <- do.call(rbind, strsplit(top_peak_factor, "_"))
  gr <- GRanges(
    seqnames = paste0("chr", parts[, 1]),
    ranges   = IRanges(start = as.numeric(parts[, 2]),
                       end   = as.numeric(parts[, 3]))
  )
  gr$peak_id <- top_peak_factor
  
  # Annotation
  peak_anno <- annotatePeak(gr,
                            TxDb         = txdb,
                            level        = "gene",
                            annoDb       = "org.Hs.eg.db",
                            tssRegion    = c(-3000, 3000),
                            verbose      = FALSE)
  
  # Conversion to dataframe 
  df <- as.data.frame(peak_anno)
  
  # reorder to have the same ordre as the one presented in the heatmap
  df <- df[match(top_peak_factor, df$peak_id), ]
  
  addWorksheet(wb, factor_name)
  writeData(wb, factor_name, df)
}

saveWorkbook(wb, file.path(outdir_graph, "peaks_anotated_by_factor.xlsx"),
             overwrite = TRUE)






