
library(mixOmics)
library(openxlsx)
source("utils.R")


input_path_rna  <- "data/RNAseq_run3_AO1_HC14_HC19_geneCounts.txt"
input_path_atac <- "data/consensus_peaks.mLb.clN.featureCounts.txt"
outdir_graph    <- "results/diablo/graph"
dir.create(outdir_graph, recursive = TRUE, showWarnings = FALSE)

#---------------------------------------------------------------------------------------------------------------------------------------
#Loading and preprocessing
#-----------------------------------------------------------------------------------------------------------------------------------------------
# ----- RNA-seq -----
rna_raw   <- read.table(input_path_rna,
                        header = TRUE, row.names = 1)
rna_raw   <- preprocess_rna(rna_raw)
rna_norm  <- normalise_rna(rna_raw)                                # VST
rna_top   <- select_top_mad(rna_norm, n_top = 10000, omic_name = "RNA")
rna       <- t(rna_top)                                       

# ----- ATAC-seq -----
atac_raw  <- read.table(input_path_atac,
                        header = TRUE, row.names = 1)

atac_raw  <- preprocess_atac(atac_raw)
atac_norm <- normalise_atac(atac_raw)                              # TMM + logCPM
atac_top  <- select_top_mad(atac_norm, n_top = 10000, omic_name = "ATAC")
atac      <- t(atac_top)

# Generate the metadata
metadata <- build_meta(t(rna))

# Order in the same way
atac <- atac[rownames(rna),]
metadata <- metadata[rownames(rna),]

#---------------------------------------------------------------------------------------------------------------------------------------
#Remove donor effect using WithinVariation 
#-----------------------------------------------------------------------------------------------------------------------------------------------
design_wv <- data.frame(sample = metadata[rownames(rna), "donor"])
rownames(design_wv) <- rownames(rna)

rna_wv  <- withinVariation(X = rna,  design = design_wv)
atac_wv <- withinVariation(X = atac, design = design_wv)

#---------------------------------------------------------------------------------------------------------------------------------------
#Diagnosis PCA
#-----------------------------------------------------------------------------------------------------------------------------------------------
plot_pca_diagnostic(t(rna),    "RNA before withinVariation")
plot_pca_diagnostic(t(atac),   "ATAC before withinVariation")
plot_pca_diagnostic(t(rna_wv), "RNA after withinVariation")
plot_pca_diagnostic(t(atac_wv),"ATAC after withinVariation")

#---------------------------------------------------------------------------------------------------------------------------------------
#Setup model mofa
#-----------------------------------------------------------------------------------------------------------------------------------------------
# Create the dataset combining rna and atacseq
X <- list(mRNA = rna_wv, atac = atac_wv)
# Create the dataset for the treatment category
Y <- as.factor(metadata[rownames(X$mRNA), "condition"])
names(Y) <- rownames(X$mRNA)
print("Number of sample per treatement :")
print(summary(Y))   

# Correlation between omics
res_pls    <- pls(X$mRNA, X$atac, ncomp = 1)
cor_blocks <- cor(res_pls$variates$X, res_pls$variates$Y)
message("Correlation between RNA-seq and Atac-seq : ", round(cor_blocks, 3))



#--------------------------------------------------------------------------------------------------------------------------------------
# Training Model
#----------------------------------------------------------------------------------------------------------------------------------------
write_metric_to_sheets <- function(wb, base_name, metric) {
  
  # If the metric is a matrix a dataframe write directly the value in the excel workshhet
  if (is.matrix(metric) || is.data.frame(metric)) {
    sheet_name <- substr(base_name, 1, 31)
    if (sheet_name %in% names(wb)) return(invisible()) # avoid duplicates
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, 
              as.data.frame(metric), 
              rowNames = TRUE,
              colNames = TRUE)
    freezePane(wb, sheet_name, firstRow = TRUE, firstCol = TRUE)
    setColWidths(wb, sheet_name, cols = 1:(ncol(metric) + 1), widths = "auto")
  }
  
  # if the metric is a list iterate over the element and write it in independant worksheet
  if (is.list(metric)) {
    for (sub_name in names(metric)) {
      new_name <- paste0(base_name, "_", sub_name)
      write_metric_to_sheets(wb, new_name, metric[[sub_name]])
    }
  }
}




value_design <- c(0.1,0.3,0.5,0.7,1)

for(value in value_design){
  design <- matrix(value, ncol = length(X), nrow = length(X), 
                   dimnames = list(names(X), names(X)))
  diag(design) <- 0
  design 
  
  diablo_explore <- block.plsda(X, Y, ncomp = 5, design = design)
  
  set.seed(123) # For reproducibility, remove for your analyses
  perf_diablo <- perf(diablo_explore,
                      validation  = 'loo',
                      dist        = c('max.dist', 'centroids.dist', 'mahalanobis.dist'),
                      progressBar = TRUE)
  
  
  #Save the plot 
  title = paste("DIABLO performance using",value,"for the design matrix")
  output_file <- paste0(outdir_graph,"/diablo_perf_",value,"_design.png")
  png(output_file, 
      width = 1200, height = 800, res = 150)
  par(oma = c(0, 0, 3, 0))
  plot(perf_diablo)
  mtext(title, side = 3, line = 0, outer = TRUE, cex = 1.3, font = 2)
  dev.off()
  
  excel_file <-  paste0(outdir_graph,"/diablo_perf_",value,"_design.xlsx")
  wb <- createWorkbook()
  error_rate_metric <- list(
    "error rate"                  = perf_diablo$error.rate, 
    "error rate per class"        = perf_diablo$error.rate.per.class,
    "Averaged Predict class"      = perf_diablo$AveragedPredict.class,
    "AveragedPredict error rate"  = perf_diablo$AveragedPredict.error.rate,
    "Weighted Predict error rate" = perf_diablo$WeightedPredict.error.rate,
    "MajorityVote error rate"     = perf_diablo$MajorityVote.error.rate,
    "WeightedVote error rate"     = perf_diablo$WeightedVote.error.rate
  )
  for (name in names(error_rate_metric)) {
    write_metric_to_sheets(wb, name, error_rate_metric[[name]])
  }
  saveWorkbook(wb, excel_file, overwrite = TRUE)

}





