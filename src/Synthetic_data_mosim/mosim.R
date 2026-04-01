library(MOSim)
library(dplyr)
library(argparse)

parser <- ArgumentParser(description = "Mosim simulation")
parser$add_argument("-a", "--donor_name", type = "character", required = TRUE,
                    help = "Name of the donor for simulation (e.g., DA1, HC7, AO1, etc.)")
parser$add_argument("-i", "--inputdir", type = "character", required = TRUE,
                    help = "Input directory for the donor data")
parser$add_argument("-o", "--outdir",   type = "character", required = TRUE,
                    help = "Output directory for saved models (.rds files)")
args <- parser$parse_args()

#-------------------------------------------------------------------------------------------------------------------------------------------
#Atac seq
#-------------------------------------------------------------------------------------------------------------------------------------------
#Paths to the data files for the ATACseq data and RNAseq data
Atacseq_featurecount_path <- paste0(args$inputdir, "/ATAC/", args$donor_name, "_peaks_consensus_peaks_featureCounts.txt")
Atacseq_annotate_path    <- paste0(args$inputdir, "/ATAC/", args$donor_name, "_peaks_consensus_peaks_annotatePeaks.txt")
rna1_path              <- paste0(args$inputdir, "/RNA/", "RNAseq_run2_HC1_HC7_geneCounts.txt")
rna2_path              <- paste0(args$inputdir, "/RNA/", "RNAseq_run3_AO1_HC14_HC19_geneCounts.txt")
#Converting the donor name to match between the ATACseq and RNAseq data (e.g., DA1 to HC1, DA2 to HC2, etc.)
if(grepl("^DA", args$donor_name)){
  individual_id <- gsub("^DA", "HC", args$donor_name)
} else {
  individual_id <- args$donor_name
}

#Creating the output directory if it does not exist
outdir_real <- paste0(args$outdir, "/real")
outdir_simu_RNA <- paste0(args$outdir, "/simu_RNA") 
outdir_simu_ATAC <- paste0(args$outdir, "/simu_ATAC")
outdir_merged <- paste0(args$outdir, "/merged_simu")

dir.create(outdir_real,      showWarnings = FALSE)
dir.create(outdir_simu_RNA,  showWarnings = FALSE)
dir.create(outdir_simu_ATAC, showWarnings = FALSE)
dir.create(outdir_merged,    showWarnings = FALSE)


#Loading the feature count file for the ATACseq data
custom_atacseq <- read.table(Atacseq_featurecount_path, header = TRUE, row.names = 1)

#Procesing the Atac data to have same column name as the RNAseq data
colnames(custom_atacseq) <- gsub("_REP.*\\.bam", "", colnames(custom_atacseq))
colnames(custom_atacseq) <- gsub("_", "", colnames(custom_atacseq))
colnames(custom_atacseq) <- gsub("DA", "HC", colnames(custom_atacseq))
colnames(custom_atacseq) <- gsub("IF$", "IFNg", colnames(custom_atacseq)) 
colnames(custom_atacseq) <- gsub("IL413IFNg", "IL413IFN", colnames(custom_atacseq))

# Create the ID of the peak by combining the Chromosome: start postion - end position
custom_atacseq <- custom_atacseq %>%
  mutate(ID = paste0(Chr, "_", Start, "_", End))
rownames(custom_atacseq) <- custom_atacseq$ID
custom_atacseq <- custom_atacseq %>% select(starts_with("HC"))

#Removing genomic patch non interval (row with names containing more the two _ since one linking the chomosome to the start and the second between start and end)
n_underscores <- stringr::str_count(rownames(custom_atacseq), "_")
valid_peaks <- rownames(custom_atacseq)[n_underscores == 2]
custom_atacseq <- custom_atacseq[valid_peaks, , drop = FALSE]

#-----------------------------------------------------------------------------------------------------------------------------------------------
#Association table
#-----------------------------------------------------------------------------------------------------------------------------------------------
annotated_atacseq <- read.table(
  Atacseq_annotate_path,
  header = TRUE,
  sep = "\t",
  quote = "",
  check.names = FALSE
)

#Processing to have the ID name with Chromosome: start postion - end position and select only the ID and gene associated
# -1 due to 1 bp offset 
#Anotation file
#Interval_152	1	11719268	11720198
#Feature count 
#Interval_152	1	11719267	11720198
#Anotation file
#Interval_15162	12	27251057	27251534
#Feature count 
#Interval_15162	12	27251056	27251534
annotated_atacseq <- annotated_atacseq %>%
  mutate(ID = paste0(Chr, "_", Start - 1, "_", End))
IDtogene <- annotated_atacseq %>% select(ID, `Entrez ID`)

colnames(IDtogene)[2] <- "Gene"
IDtogene <- IDtogene %>%
  arrange(match(ID, rownames(custom_atacseq)))
intersect(rownames(custom_atacseq), IDtogene$ID) %>% length()
IDtogene <- IDtogene %>% filter(ID %in% valid_peaks)

#-----------------------------------------------------------------------------------------------------------------------------------------------
#RNAseq
#-----------------------------------------------------------------------------------------------------------------------------------------------
rna1 = read.table(rna1_path, header = TRUE, row.names = 1)
rna2 = read.table(rna2_path, header = TRUE, row.names = 1)

#Verify if same order row before merging
stopifnot(rownames(rna1) == rownames(rna2)) 
custom_rnaseq <- cbind(rna1, rna2)

#Delete the unmapp, multimap row
custom_rnaseq <- custom_rnaseq[!grepl("^N_", rownames(custom_rnaseq)), ]

#-------------------------------------------------------------------------------------------------------------------------------------------------
#Simulation
#-------------------------------------------------------------------------------------------------------------------------------------------------

cat("List of treatment for the ATACseq", sort(colnames(custom_atacseq)), "\n")
column_donor <- grep(paste0("^", individual_id), colnames(custom_rnaseq), value = TRUE)
cat("List of treatment for the RNAseq", sort(column_donor), "\n")


#Save the real data for the donor
write.table(custom_atacseq, paste0(outdir_real, "/", individual_id, "_real_ATAC.tsv"), sep = "\t", quote = FALSE, col.names = NA)
write.table(custom_rnaseq[, sort(column_donor)], paste0(outdir_real, "/", individual_id, "_real_RNA.tsv"), sep = "\t", quote = FALSE, col.names = NA)

for(current_ind_treat in colnames(custom_atacseq)){
  if(!current_ind_treat %in% colnames(custom_rnaseq)){
    warning(paste("No:", current_ind_treat, "treatment found in the RNAseq data"))
    next
  }
  rna_seed <- custom_rnaseq %>% select(all_of(current_ind_treat))
  colnames(rna_seed) <- "Counts"
  
  atac_seed <- custom_atacseq %>% select(all_of(current_ind_treat))
  colnames(atac_seed) <- "Counts"
  
  # Library size in millions for depth parameter
  rna_depth  <- sum(rna_seed$Counts)  / 1e6
  atac_depth <- sum(atac_seed$Counts) / 1e6
  cat("Library size RNA:", rna_depth, "M | ATAC:", atac_depth, "M for", current_ind_treat, "\n")
  
  #Creation of the Mosim Data format for the simulation
  rna_data <- omicData("RNA-seq", data = rna_seed)
  
  atac_data <- omicData(
    "DNase-seq",
    data            = atac_seed,
    associationList = IDtogene)
  
  # Simulation
  simulation <- mosim(
    omics        = c("RNA-seq", "DNase-seq"),
    times        = 0,
    numberGroups = 2,
    numberReps   = 10,
    diffGenes    = 0.0001,
    omicsOptions = list(
      "RNA-seq" = list(
        data            = rna_seed,
        depth           = rna_depth
      ),
      "DNase-seq" = list(
        data            = atac_seed,
        idToGene        = IDtogene,
        depth           = atac_depth
      )
    )
  )
  
  print(simulation@simulators$SimRNAseq@depth)
  print(simulation@simulators$SimDNaseseq@depth)
  # 4. Verification correct number and name genes and group and setting of the simulation
  print(omicSettings(simulation))
  
  experimentalDesign(simulation)
  
  # 5. Retrieve the simulated data and save them
  dataRNAseq  <- omicResults(simulation, "RNA-seq")
  dataATACseq <- omicResults(simulation, "DNase-seq")
  
  #Selecting only group 1 and renaming the column with the individual and treatment name
  dataRNAseq  <- dataRNAseq  %>% select(!starts_with("Group2"))
  dataATACseq <- dataATACseq %>% select(!starts_with("Group2"))
  colnames(dataRNAseq)  <- gsub("Group1.Time0", current_ind_treat, colnames(dataRNAseq))
  colnames(dataATACseq) <- gsub("Group1.Time0", current_ind_treat, colnames(dataATACseq))
  
  # Saving the files
  write.table(dataRNAseq,  paste0(outdir_simu_RNA, "/", current_ind_treat, "RNA.tsv"),  sep = "\t", quote = FALSE, col.names = NA)
  write.table(dataATACseq, paste0(outdir_simu_ATAC, "/", current_ind_treat, "ATAC.tsv"), sep = "\t", quote = FALSE, col.names = NA)
  cat("Simulation of", current_ind_treat, "complete\n")
}

#-------------------------------------------------------------------------------------------------------------------------------------------------
#Merge simulated files
#-------------------------------------------------------------------------------------------------------------------------------------------------
rna_files <- list.files(outdir_simu_RNA,
                        pattern = paste0("^", individual_id, ".*RNA\\.tsv$"),
                        full.names = TRUE)
cat("RNA files found:\n")
cat(basename(rna_files), sep = "\n")

sim_rna <- do.call(cbind, lapply(rna_files, function(f) {
  read.table(f, header = TRUE, row.names = 1, sep = "\t")
}))
cat("RNA merged:", nrow(sim_rna), "x", ncol(sim_rna), "\n")

atac_files <- list.files(outdir_simu_ATAC,
                         pattern = paste0("^", individual_id, ".*ATAC\\.tsv$"),
                         full.names = TRUE)
cat("ATAC files found:\n")
cat(basename(atac_files), sep = "\n")

sim_atac <- do.call(cbind, lapply(atac_files, function(f) {
  read.table(f, header = TRUE, row.names = 1, sep = "\t")
}))
cat("ATAC merged:", nrow(sim_atac), "x", ncol(sim_atac), "\n")

write.table(sim_rna,
            paste0(outdir_merged, "/", individual_id, "_RNA_merged.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
write.table(sim_atac,
            paste0(outdir_merged, "/", individual_id, "_ATAC_merged.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
cat("Merge complete for", individual_id, "\n")