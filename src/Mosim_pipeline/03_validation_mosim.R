#' Validation of mosim simulated data using countsimQC
#'
#' This script takes as input the simulated of one donor and the real data of the corresponding donor, for both RNA and ATAC. 
#' It then creates DESeqDataSet objects for the simulated and real data, runs countsimQC to compare the distributions of counts, and generates an HTML report.
#' @param donor_ID ID of the donor for simulation in the format HC followed by a number (e.g., HC7, HC19, etc.)
#' @param inputdir Input directory containing the real and simulated data for the specified donor. The script expects the following files in this directory:
#'   - real/HC19_real_RNA.tsv
#'   - merged_simu/HC19_RNA_merged.tsv
#'   - real/HC19_real_ATAC.tsv
#'   - merged_simu/HC19_ATAC_merged.tsv
#' @param outdir Output directory where the countsimQC HTML reports will be saved. The reports will be named as follows:
#'   - HC19_1rep_vs_allreal_RNA.html
#'   - HC19_1rep_vs_allreal_ATAC.html

#' @examples
#' Rscript validation_mosim.R \
#'    --donor_ID "HC19" \
#'    --inputdir "/path/to/simulated_data" \
#'    --outdir "/path/to/validation_results"

library(DESeq2)
library(countsimQC)
library(argparse)

parser <- ArgumentParser(description = "Validation of simulated data using countsimQC")
parser$add_argument("-a", "--donor_ID", type = "character", required = TRUE,
                    help = "ID of the donor for simulation (e.g., HC1, HC7, HC19, etc.)")
parser$add_argument("-i", "--inputdir", type = "character", required = TRUE,
                    help = "Input directory for the donor data")
parser$add_argument("-o", "--outdir", type = "character", required = TRUE,
                    help = "Output directory for saved models (.rds files)")
parser$add_argument("-r", "--replicate", type = "integer", default = 1,
                    help = "Replicate number to compare (default: 1)")
args <- parser$parse_args()


individual_id <- args$donor_ID
data_dir      <- args$inputdir
output_dir    <- args$outdir
replicate_number <- args$replicate

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


load_data <- function(modality) {
  #' Function to load the real and simulated data for a given modality (RNA or ATAC) for the specified donor.
  #' @param modality Modality to load ("RNA" or "ATAC")
  #' @return A list containing the real and simulated data frames for the specified modality.
  #' @examples
  #' data_lists <- load_data("RNA")
  #' real_rna <- data_lists$real
  #' sim_rna <- data_lists$simulated
  real_path <- file.path(data_dir, "real", paste0(individual_id, "_real_", modality, ".tsv"))
  sim_path  <- file.path(data_dir, "merged_simu", paste0(individual_id, "_", modality, "_merged.tsv"))

  if (!file.exists(real_path)) {
    stop("Real data file not found: ", real_path)
  }
  if (!file.exists(sim_path)) {
    stop("Simulated data file not found: ", sim_path)
  }

  real_data <- read.table(real_path, header = TRUE, row.names = 1, sep = "\t")
  sim_data  <- read.table(sim_path,  header = TRUE, row.names = 1, sep = "\t")
  cat("Dimensions of the data for", modality, ":\n")
  cat("Real", modality, ":", nrow(real_data), "x", ncol(real_data), "\n")
  cat("Simulated", modality, ":", nrow(sim_data), "x", ncol(sim_data), "\n")
  return(list(real = real_data, simulated = sim_data))
}

retrieve_replicate <- function(sim_data, modality, replicate = 1) {
  #' Function to retrieve the specified replicate from the simulated data for a given modality.
  #' @param sim_data Data frame of the simulated data for the specified modality.
  #' @param modality Modality to retrieve ("RNA" or "ATAC").
  #' @param replicate Replicate number to retrieve.
  #' @return A data frame containing the specified replicate for the given modality.

  replicate_suffix <- paste0(".Rep", replicate, "$")
  rep_cols <- grep(replicate_suffix, colnames(sim_data), value = TRUE)

  if (length(rep_cols) == 0) {
    stop("No columns matching pattern '", replicate_suffix, "' found in ", modality,
         ". Available columns: ", paste(head(colnames(sim_data)), collapse = ", "))
  }

  sim_rep <- sim_data[, rep_cols, drop = FALSE]
  colnames(sim_rep) <- gsub(replicate_suffix, "", colnames(sim_rep))
  cat("Extracted Rep", replicate, "columns for", modality, ":", length(rep_cols), "\n")
  return(sim_rep)
}

generate_deseq_dataset <- function(count_data, modality) {
  #' Function to create a DESeqDataSet object from the given count data for a specified modality.
  #' @param count_data Data frame of count data (features x samples) for the specified modality.
  #' @param modality Modality for which the DESeqDataSet is being created ("RNA" or "ATAC").
  #' @return A DESeqDataSet object created from the count data.
  dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(count_data)),
    colData   = data.frame(condition = factor(colnames(count_data)),
                           row.names = colnames(count_data)),
    design    = ~1
  )
  cat("Created DESeqDataSet for", modality, "with dimensions:", dim(dds), "\n")
  return(dds)
}
generate_countsimqc_report <- function(dds_real, dds_sim, modality, output_dir, individual_id,replicate = 1) {
  #' Function to run countsimQC and generate an HTML report comparing the real and simulated data for a given modality.
  #' @param dds_real DESeqDataSet object for the real data.
  #' @param dds_sim DESeqDataSet object for the simulated data.
  #' Parameter for saving the report filename:
  #' @param modality Modality being compared ("RNA" or "ATAC").
  #' @param output_dir Directory where the countsimQC report will be saved.
  #' @param individual_id ID of the donor being analyzed.
  #' @param replicate Replicate number being compared.

  output_file <- paste0(individual_id,"_",replicate,"_vs_allreal_", modality, ".html")
  countsimQCReport(
    ddsList             = list(Real = dds_real, Simulated = dds_sim),
    outputFile          = output_file,
    outputDir           = output_dir,
    outputFormat        = "html_document",
    showCode            = FALSE,
    calculateStatistics = TRUE
  )
  cat("Generated countsimQC report for", modality, ":", file.path(output_dir, output_file), "\n")
}
run_validation <- function(modality) {
  #' Function to run the full validation process for a given modality, including loading data, retrieving the specified replicate,
  #' creating DESeqDataSet objects, and generating the countsimQC report.
  #' @param modality Modality to validate ("RNA" or "ATAC").
  cat("\n=== Processing", modality, "===\n")
  data_lists <- load_data(modality)
  sim_rep <- retrieve_replicate(data_lists$simulated, modality, replicate_number)
  dds_real <- generate_deseq_dataset(data_lists$real, modality)
  dds_sim <- generate_deseq_dataset(sim_rep, paste("Simulated", modality))
  generate_countsimqc_report(dds_real, dds_sim, modality, output_dir, individual_id, replicate_number)
}

# Run the validation for both RNA and ATAC
run_validation("RNA")
run_validation("ATAC")

cat("Validation completed for donor", individual_id, "with replicate", replicate_number, "\n")