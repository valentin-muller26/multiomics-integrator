library(DESeq2)
library(countsimQC)
library(argparse)

parser <- ArgumentParser(description = "Validation of simulated data using countsimQC")
parser$add_argument("-a", "--donor_name", type = "character", required = TRUE,
                    help = "Name of the donor for simulation (e.g., DA1, HC7, AO1, etc.)")
parser$add_argument("-i", "--inputdir", type = "character", required = TRUE,
                    help = "Input directory for the donor data")
parser$add_argument("-o", "--outdir", type = "character", required = TRUE,
                    help = "Output directory for saved models (.rds files)")
args <- parser$parse_args()

individual_id <- args$donor_name
data_dir      <- args$inputdir
output_dir    <- args$outdir

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Chargement des données
real_rna <- read.table(paste0(data_dir, "/real/", individual_id, "_real_RNA.tsv"),
                       header = TRUE, row.names = 1, sep = "\t")
sim_rna  <- read.table(paste0(data_dir, "/merged_simu/", individual_id, "_RNA_merged.tsv"),
                       header = TRUE, row.names = 1, sep = "\t")

cat("Real RNA :", nrow(real_rna), "x", ncol(real_rna), "\n")
cat("Sim RNA  :", nrow(sim_rna),  "x", ncol(sim_rna),  "\n")



#-------------------------------------------------------------------------------------------------------------------------------------------------
# Comparaison 2 : tous les simulés (80 colonnes) vs tous les réels (8 colonnes)
#-------------------------------------------------------------------------------------------------------------------------------------------------
countsimQCReport(
  ddsList             = list(Real = real_rna, Simulated = sim_rna),
  outputFile          = paste0(individual_id, "_allsimu_vs_allreal_countsimQC.html"),
  outputDir           = output_dir,
  outputFormat        = "html_document",
  showCode            = FALSE,
  calculateStatistics = FALSE
)
cat("Comparaison 2 done: tous simulés vs tous réels\n")

#-------------------------------------------------------------------------------------------------------------------------------------------------
# Comparaison 3 : 1 réplicat simulé par traitement (8 colonnes) vs réels (8 colonnes)
#-------------------------------------------------------------------------------------------------------------------------------------------------
rep1_cols <- grep("\\.Rep1$", colnames(sim_rna), value = TRUE)
sim_rep1  <- sim_rna[, rep1_cols, drop = FALSE]
colnames(sim_rep1) <- gsub("\\.Rep1$", "", colnames(sim_rep1))

dds_real <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(real_rna)),
  colData   = data.frame(condition = factor(colnames(real_rna)),
                         row.names = colnames(real_rna)),
  design    = ~1
)

dds_sim <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(sim_rep1)),
  colData   = data.frame(condition = factor(colnames(sim_rep1)),
                         row.names = colnames(sim_rep1)),
  design    = ~1
)

countsimQCReport(
  ddsList             = list(Real = dds_real, Simulated = dds_sim),
  outputFile          = paste0(individual_id, "_1rep_vs_allreal_countsimQC.html"),
  outputDir           = output_dir,
  outputFormat        = "html_document",
  showCode            = FALSE,
  calculateStatistics = FALSE
)
cat("Comparaison 3 done: 1 replicat par traitement vs réels\n")