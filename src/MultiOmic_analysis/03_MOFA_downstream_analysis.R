#' Downstream analysis of a trained MOFA2 model (RNA-seq + ATAC-seq)
#' @description 
#' This script performs post-training analysis of a MOFA2 model
#   1. Load the trained model and pre-computed RNA / ATAC data
#   2. Visualise explained variance and correlations with covariates
#   3. Explore latent factors (1D and 2D) by condition / donor
#   4. Plot heatmaps of top genes / peaks for condition-associated factors
#   5. Run gene set enrichment analysis (GSEA) on the RNA view
#   6. Inspect specific genes of interest in the factor space
#   7. Annotate ATAC-seq peaks with genomic features via ChIPseeker
#'@param outdir_model Output directory for the trained model and preprocessed matrices (.rds)
#'@param outdir_graph Output directory for the diagnostic plots (PDF)
 

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
library(argparse)

if ("package:data.table" %in% search()) {
  detach("package:data.table", unload = TRUE, character.only = TRUE)
}

parser <- ArgumentParser(description = " Downstream analysis of a trained MOFA2 model (RNA-seq + ATAC-seq)")
parser$add_argument("-f", "--nb_factor",
                    type = "integer", default = 5,
                    help = "Maximum number of latent factors to fit. (default: 5)")

parser$add_argument("-om", "--outdir_model",
                    type = "character", required = TRUE,
                    help = "Output directory where the trained MOFA model (.hdf5) will be saved.")

parser$add_argument("-og", "--outdir_graph",
                    type = "character", required = TRUE,
                    help = "Output directory for plots generated during the training (elbo and comparison models).")
args <- parser$parse_args()


# Source the functions for MOSim simulation and data processing
path_scripts <- Sys.getenv("PATH_SCRIPTS")
if (path_scripts == "") {
    stop("Environment variable PATH_SCRIPTS is not set. Source 00_config.sh before running.")
}
source(file.path(path_scripts, "utils.R"))
source(file.path(path_scripts, "MOFA_function.R"))


# Set the output directory
outdir_model    <- args$outdir_model
outdir_graph    <-  args$outdir_graph
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
sig_factors_names <-  c("Factor1", "Factor3","Factor5")
list_factor_condition <- as.numeric(gsub("Factor", "", sig_factors_names))
genes_of_interest <- list("ENSG00000232810","ENSG00000136634") #list of the gene of interest for the Specific gene anaylsis





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

# Anaylsis using the ENSEMBL_ID
top_feature <- heatmap_analysis(MOFAobject.trained,list_factor_condition)
top_peaks <-top_feature[["top_peaks"]]
top_gene <- top_feature[["top_gene"]]

#Analysis using the gene_name
heatmap_analysis(MOFAobject.trained,list_factor_condition,gene_name=TRUE)


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
#    reactomeGS : Reactome pathway database
#    MSigDB C2  : curated gene sets
#    MSigDB C5  : Gene Ontology gene sets
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
top_peaks
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



