

library(MOSim)
library(dplyr)
library(argparse)

parser <- ArgumentParser(description = "Mosim simulation")
parser$add_argument("-a", "--donor_ID", type = "character", required = TRUE,
                    help = "Name of the donor for simulation (e.g., DA1, HC7, AO1, etc.)")
parser$add_argument("-i", "--inputdir", type = "character", required = TRUE,
                    help = "Input directory for the donor data")
parser$add_argument("-o", "--outdir",   type = "character", required = TRUE,
                    help = "Output directory for saved models (.rds files)")
args <- parser$parse_args()

# Set seed for reproducibility
set.seed(42)

# Source the functions for MOSim simulation and data processing
path_scripts <- Sys.getenv("PATH_SCRIPTS")
if (path_scripts == "") {
    stop("Environment variable PATH_SCRIPTS is not set. Source 00_config.sh before running.")
}
source(file.path(path_scripts, "mosim_functions.R"))


# Path to the input files for the specified donor
Atacseq_featurecount_path <- file.path(args$inputdir, "ATAC", paste0(args$donor_ID, "_peaks_consensus_peaks_featureCounts.txt"))
Atacseq_annotate_path     <- file.path(args$inputdir, "ATAC", paste0(args$donor_ID, "_peaks_consensus_peaks_annotatePeaks.txt"))
rna1_path                 <- file.path(args$inputdir, "RNA", "RNAseq_run2_HC1_HC7_geneCounts.txt")
rna2_path                 <- file.path(args$inputdir, "RNA", "RNAseq_run3_AO1_HC14_HC19_geneCounts.txt")

# Path the output directories for the simulated data 
outdir_real      <- file.path(args$outdir, "real")
outdir_simu_RNA  <- file.path(args$outdir, "simu_RNA")
outdir_simu_ATAC <- file.path(args$outdir, "simu_ATAC")
outdir_merged    <- file.path(args$outdir, "merged_simu")

# Create output directories if they don't exist
dir.create(outdir_real,      showWarnings = FALSE, recursive = TRUE)
dir.create(outdir_simu_RNA,  showWarnings = FALSE, recursive = TRUE)
dir.create(outdir_simu_ATAC, showWarnings = FALSE, recursive = TRUE)
dir.create(outdir_merged,    showWarnings = FALSE, recursive = TRUE)

# Ensure that the donor ID is consistent between RNA and ATAC data (replace "DA" with "HC" if needed)
if(grepl("^DA", args$donor_ID)){
  individual_id <- gsub("^DA", "HC", args$donor_ID)
} else {
  individual_id <- args$donor_ID
}

#--------------------------------------------------------------------------------------------------------------------------------------------------
# Load and preprocess the different omics and annotation files for the specified donor
#--------------------------------------------------------------------------------------------------------------------------------------------------
# Load and preprocess the ATAC-seq data
custom_atacseq <- read.table(Atacseq_featurecount_path, header = TRUE, row.names = 1)
atac_result <- preprocess_ATACseq(custom_atacseq)
custom_atacseq <- atac_result$data
valid_peaks    <- atac_result$valid_peaks

# load the annotation file and create a mapping of peak IDs to gene IDs for the ATAC-seq data
annotated_atacseq <- read.table(
  Atacseq_annotate_path,
  header = TRUE,
  sep = "\t",
  quote = "",
  check.names = FALSE
)
IDtogene <- preprocess_annotation(annotated_atacseq, valid_peaks, custom_atacseq)
 
# Load and merge the RNA-seq data for the specified donor
custom_rnaseq <- load_merge_rna(rna1_path, rna2_path)


#-------------------------------------------------------------------------------------------------------------------------------------------------
# Preprocess the real data (filtering lowly expressed genes and lowly accessible peaks) and save the filtered real data for the specified donor
#-------------------------------------------------------------------------------------------------------------------------------------------------
cat("List of treatment for the ATACseq", sort(colnames(custom_atacseq)), "\n")
column_donor <- grep(paste0("^", individual_id), colnames(custom_rnaseq), value = TRUE)
cat("List of treatment for the RNAseq", sort(column_donor), "\n")

cat("Number of peaks in ATACseq before filtering:", nrow(custom_atacseq), "\n")
cat("Number of genes in RNAseq before filtering:", nrow(custom_rnaseq), "\n")

# Filter out lowly expressed genes and lowly accessible peaks based on mean and variance thresholds 
custom_rnaseq <-filter_low_count(custom_rnaseq, mean_threshold = 0, var_threshold = 0.1)
custom_atacseq <- filter_low_count(custom_atacseq, mean_threshold = 0, var_threshold = 0.1)
cat("Number of peaks in ATACseq after filtering:", nrow(custom_atacseq), "\n")
cat("Number of genes in RNAseq after filtering:", nrow(custom_rnaseq), "\n")  

#Save the real data for the donor
write.table(custom_atacseq,
            file.path(outdir_real, paste0(individual_id, "_real_ATAC.tsv")),
            sep = "\t", quote = FALSE, col.names = NA)
write.table(custom_rnaseq[, sort(column_donor)],
            file.path(outdir_real, paste0(individual_id, "_real_RNA.tsv")),
            sep = "\t", quote = FALSE, col.names = NA)

#-------------------------------------------------------------------------------------------------------------------------------------------------
# Simulate one replicate per treatment (8 columns) 
#-------------------------------------------------------------------------------------------------------------------------------------------------
simulation_one_donor(
    rna_df           = custom_rnaseq,
    atac_df          = custom_atacseq,
    IDtogene         = IDtogene,
    outdir_simu_RNA  = outdir_simu_RNA,
    outdir_simu_ATAC = outdir_simu_ATAC
)

#-------------------------------------------------------------------------------------------------------------------------------------------------
#Merge simulated files
#-------------------------------------------------------------------------------------------------------------------------------------------------
merge_simulated_files(outdir_simu_RNA,  individual_id, "RNA",  outdir_merged)
merge_simulated_files(outdir_simu_ATAC, individual_id, "ATAC", outdir_merged)
cat("Merge complete for", individual_id, "\n")

#-------------------------------------------------------------------------------------------------------------------------------------------------
# Cleanup: remove per-treatment intermediate files (only keep the merged ones)
#-------------------------------------------------------------------------------------------------------------------------------------------------
unlink(outdir_simu_RNA,  recursive = TRUE)
unlink(outdir_simu_ATAC, recursive = TRUE)
cat("Removed intermediate directories:", outdir_simu_RNA, "and", outdir_simu_ATAC, "\n")