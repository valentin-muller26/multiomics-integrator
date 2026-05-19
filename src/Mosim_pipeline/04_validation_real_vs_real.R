#' Validation of real data using countsimQC
#' This script compares the real data of two real donors for a given modality (RNA or ATAC) using countsimQC. 
#' It loads the real data for both donors, compares the number of features and their overlap.
#' Then, it creates DESeqDataSet objects for both donors and generates a countsimQC report comparing the two datasets.
#' This script allow to compare real vs real data to have a baseline of the variability between real donors,
#' which can be used to interpret the results of the comparison between simulated and real data.
#' @param input_path Path to the input data directory containing the real data files (e.g., /data/users/vmuller/0_master_thesis/data/data_interleukines/simulated_data/real)
#' @param donor_1 Name of the first donor for comparison (e.g., HC17, HC19, etc.)
#' @param donor_2 Name of the second donor for comparison (e.g., HC17, HC19, etc.)
#' @param modality Modality to compare (RNA or ATAC)
#' @param outdir Output directory for the countsimQC report (.html file)
#' @examples
#' Rscript validation_reel_atac.R \
#'    --input_path "/data/users/vmuller/0_master_thesis/data/data_interleukines/simulated_data/real" \
#'    --donor_1 "HC17" \
#'    --donor_2 "HC19" \
#'    --modality "ATAC" \
#'    --outdir "/data/users/vmuller/0_master_thesis/data/data_interleukines/real_validation_results" 

library(DESeq2)
library(countsimQC)
library(argparse)

parser <- ArgumentParser(description = "Validation of real data using countsimQC")
parser$add_argument("-i", "--input_path", type = "character", required = TRUE,
                    help = "Path to the input data directory containing the real data files (e.g., /data/users/vmuller/0_master_thesis/data/data_interleukines/simulated_data/real)")
parser$add_argument("-a", "--donor_1", type = "character", required = TRUE,
                    help = "Name of the first donor for comparison (e.g., DA1, HC7, AO1, etc.)")
parser$add_argument("-b", "--donor_2", type = "character", required = TRUE,
                    help = "Name of the second donor for comparison (e.g., DA2, HC8, AO2, etc.)")
parser$add_argument("-m", "--modality", type = "character", required = TRUE,
                    help = "Modality to compare (RNA or ATAC)")
parser$add_argument("-o", "--outdir",   type = "character", required = TRUE,
                    help = "Output directory for the countsimQC report  (.html file)")
args <- parser$parse_args()


load_data <- function(donor_name, modality) {
  #' Function to load the real data for a given donor and modality (RNA or ATAC).
  #' @param donor_name Name of the donor (e.g., HC17, HC19, etc.)
  #' @param modality Modality to load ("RNA" or "ATAC")
  #' @return A data frame containing the real data for the specified donor and modality.
  #' @examples
  #' real_data <- load_data("HC17", "RNA")
  data_path <- file.path(args$input_path,
                         paste0(donor_name, "_real_", modality, ".tsv"))
  if (!file.exists(data_path)) {
    stop("Data file not found: ", data_path)
  }
  data <- read.table(data_path, header = TRUE, row.names = 1, sep = "\t")
  cat("Dimensions of the real data for", donor_name, "and", modality, ":", nrow(data), "x", ncol(data), "\n")
  return(data)
}

comparison_donors <-  function(real_data_1, real_data_2, donor_1, donor_2, modality) {
  #' Function to compare the real data of two donors for the number of features (genes or peaks) and the overlap between them.
  #' @param donor_1 Name of the first donor (e.g., HC17, HC19, etc.)
  #' @param donor_2 Name of the second donor (e.g., HC17, HC19, etc.)
  #' @param modality Modality to compare ("RNA" or "ATAC")
  #' @examples
  #' comparison_donors("HC17", "HC19", "ATAC")

  # Peaks in each dataset
  feature_donor_1 <- rownames(real_data_1)
  feature_donor_2 <- rownames(real_data_2)


  # Common peaks
  feature_common <- intersect(feature_donor_1, feature_donor_2)

  # Unique peaks
  feature_only_donor_1 <- setdiff(feature_donor_1, feature_donor_2)
  feature_only_donor_2  <- setdiff(feature_donor_2, feature_donor_1)

  # Summary
  cat("Comparison of", modality, "features between", donor_1, "and", donor_2, ":\n")
  feature_name <- ifelse(modality == "RNA", "Genes", "Peaks")
  cat (feature_name, donor_1, ":", length(feature_donor_1), "\n")
  cat (feature_name, donor_2, ":", length(feature_donor_2), "\n")
  cat("Common", feature_name, ":", length(feature_common), "\n")
  cat("Unique", feature_name, "in", donor_1, ":", length(feature_only_donor_1), "\n")
  cat("Unique", feature_name, "in", donor_2, ":", length(feature_only_donor_2), "\n")
  cat("% common", feature_name, "out of", donor_1, ":", round(length(feature_common) / length(feature_donor_1) * 100, 2), "%\n")
  cat("% common", feature_name, "out of", donor_2, ":", round(length(feature_common) / length(feature_donor_2) * 100, 2), "%\n")
} 

create_DESeqDataSet <- function(data) {
  #' Function to create a DESeqDataSet object from the real data for a given donor
  #' @param data Data frame of the real data for the specified donor and modality.
  #' @return A DESeqDataSet object containing the count data and colData for the specified donor and modality.
  #' @examples
  #' dds <- create_DESeqDataSet(real_data)
  dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(data)),
    colData   = data.frame(condition = factor(colnames(data)),
                           row.names = colnames(data)),
    design    = ~1
  )
  return(dds)
}

generate_report_countsimQC <- function(dds_real_1, dds_real_2, donor_1, donor_2, modality) {
  #' Function to generate a countsimQC report comparing the real data of two donors for a given modality.
  #' @param dds_real_1 DESeqDataSet object for the first donor's real data.
  #' @param dds_real_2 DESeqDataSet object for the second donor's real data.
  #' @param donor_1 Name of the first donor (e.g., HC17, HC19, etc.)
  #' @param donor_2 Name of the second donor (e.g., HC17, HC19, etc.)
  #' @param modality Modality to compare ("RNA" or "ATAC")
  #' @examples
  #' generate_report_countsimQC(dds_real_1, dds_real_2, "HC17", "HC19", "ATAC")

  report_file <- paste0(donor_1, "_vs_", donor_2, "_real_", modality, ".html")
  countsimQCReport(
    ddsList             = list(Real_Donor1 = dds_real_1, Real_Donor2 = dds_real_2),
    outputFile          = report_file,
    outputDir           = args$outdir,
    outputFormat        = "html_document",
    showCode            = FALSE,
    calculateStatistics = TRUE
  )
  cat("Comparison", donor_1, "vs", donor_2, "for", modality, "done. Report saved to:", file.path(args$outdir, report_file), "\n")
}

#-------------------------------------------------------------------------------
# Run the validation steps
#-------------------------------------------------------------------------------
donor_1  <- args$donor_1
donor_2  <- args$donor_2
modality <- args$modality

# Create output directory if it doesn't exist
dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)

real_data_1 <- load_data(donor_1, modality)
real_data_2 <- load_data(donor_2, modality)
comparison_donors(real_data_1, real_data_2, donor_1, donor_2, modality)
dds_1 <- create_DESeqDataSet(real_data_1)
dds_2 <- create_DESeqDataSet(real_data_2)
generate_report_countsimQC(dds_1, dds_2, donor_1, donor_2, modality)