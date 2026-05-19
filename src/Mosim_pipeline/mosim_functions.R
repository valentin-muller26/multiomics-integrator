preprocess_ATACseq <- function(atac_df) {
    #' Function to preprocess the ATAC-seq data by cleaning column names, creating unique peak IDs composed of chromosome, start, and end coordinates, and filtering valid peaks.
    #' This function align the name of the samples (in the column names) in the ATAC-seq data with the RNA-seq data and select only the column with the ID and feature count    
    #' @param custom_atacseq Data frame of the raw ATAC-seq feature counts with columns for chromosome, start, end, and sample counts.
    #' @return A preprocessed data frame of ATAC-seq counts with unique peak IDs as row names and cleaned sample names as column names.
    #' @examples
    #' preprocessed_atacseq <- preprocess_ATACseq(raw_atacseq)

    # Clean column names to match RNA-seq sample names
    colnames(atac_df) <- gsub("_REP.*\\.bam", "", colnames(atac_df))
    colnames(atac_df) <- gsub("_", "", colnames(atac_df))
    colnames(atac_df) <- gsub("DA", "HC", colnames(atac_df))
    colnames(atac_df) <- gsub("IF$", "IFNg", colnames(atac_df)) 
    colnames(atac_df) <- gsub("IL413IFNg", "IL413IFN", colnames(atac_df))

    # Create unique peak IDs by combining chromosome, start, and end coordinates
    custom_atacseq <- atac_df %>%
      mutate(ID = paste0(Chr, "_", Start, "_", End))
    rownames(custom_atacseq) <- custom_atacseq$ID
    custom_atacseq <- custom_atacseq %>% select(starts_with("HC"))

    # Filter valid peaks (those with exactly 2 underscores in their ID)
    n_underscores <- stringr::str_count(rownames(custom_atacseq), "_")
    valid_peaks <- rownames(custom_atacseq)[n_underscores == 2]
    custom_atacseq <- custom_atacseq[valid_peaks, , drop = FALSE]
    return(list(data = custom_atacseq, valid_peaks = valid_peaks))
}

preprocess_annotation <- function(annot_df, valid_peaks, custom_atacseq) {
  #' @param annot_df Raw ATAC-seq annotation data frame.
  #' @param valid_peaks Vector of valid peak IDs (output from preprocess_ATACseq).
  #' @param custom_atacseq Preprocessed ATAC-seq data frame.
  #' @return Named data frame mapping peak IDs to Entrez gene IDs.

  annot_df <- annot_df %>%
    mutate(ID = paste0(Chr, "_", Start - 1, "_", End))

  IDtogene <- annot_df %>% select(ID, `Entrez ID`)
  colnames(IDtogene)[2] <- "Gene"

  IDtogene <- IDtogene %>%
    arrange(match(ID, rownames(custom_atacseq))) %>%
    filter(ID %in% valid_peaks)

  return(IDtogene)
}

load_merge_rna <- function(rna1_path, rna2_path){
    #' Function to load and merge RNA-seq data from two runs, ensuring that the row names (genes) match and filtering out any unwanted features (eg not mapped genes).
    #' @param rna1_path File path to the first RNA-seq data file (e.g., "RNAseq_run2_HC1_HC7_geneCounts.txt").
    #' @param rna2_path File path to the second RNA-seq data file (e.g., "RNAseq_run3_AO1_HC14_HC19_geneCounts.txt").
    #' @return A merged data frame of RNA-seq counts with genes as row names and samples as column names, filtered to remove any features that start with "N_" (not mapped genes).
    #' @examples
    #' merged_rna <- load_merge_rna(rna1_path, rna2_path)

    rna1 = read.table(rna1_path, header = TRUE, row.names = 1)
    rna2 = read.table(rna2_path, header = TRUE, row.names = 1)

    stopifnot(identical(rownames(rna1), rownames(rna2)))
    custom_rnaseq <- cbind(rna1, rna2)
    custom_rnaseq <- custom_rnaseq[!grepl("^N_", rownames(custom_rnaseq)), ]
    return(custom_rnaseq)
}

filter_low_count <- function(omics_df, mean_threshold = 0, var_threshold = 0.1) {
    #' Function to filter out lowly expressed genes or lowly accessible peaks based on mean and variance thresholds.
    #' @param omics_df Data frame of count data (features x samples) for either RNA-seq or ATAC-seq.
    #' @param mean_threshold Minimum mean count across samples for a feature to be retained (default is 0).
    #' @param var_threshold Minimum variance of log-transformed counts across samples for a feature to be retained (default is 0.1).
    #' @return A filtered data frame of count data with lowly expressed genes or lowly accessible peaks removed based on the specified thresholds.
    #' @examples
    #' filtered_df <- filter_low_count(omics_df, mean_threshold = 0, var_threshold = 0.1)
    omics_mask <- rowMeans(omics_df) > mean_threshold  &
            apply(log1p(omics_df), 1, var) >= var_threshold
    omics_df <- omics_df[omics_mask, ]
    return(omics_df)
}

prepare_seed_data <- function(df, treatment) {
    #' Extract one treatment column from a count data frame and prepare it as MOSim seed input.
    #' MOSim expects a single-column data frame with the column named "Counts".
    #' @param df Data frame of count data (features x samples).
    #' @param treatment Name of the column to extract.
    #' @return A list with two elements:
    #'   - data: single-column data frame named "Counts".
    #'   - depth: library depth in millions of reads.

    seed <- df %>% select(all_of(treatment))
    colnames(seed) <- "Counts"
    depth <- sum(seed$Counts) / 1e6
    return(list(data = seed, depth = depth))
}

postprocess_simulated_data <- function(simulation, omic_name, treatment) {
    #' Function to post-process the simulated data from MOSim by selecting only the columns corresponding to Group1 (treatment without the logfold change applied)
    #' and renaming them to match the original treatment name.
    #' @param simulation The output object from the MOSim simulation.
    #' @param omic_name Name of the omic modality to extract ("RNA-seq" or "DNase-seq").
    #' @param treatment Name of the original treatment to use for renaming the columns in the output.
    #' @return A data frame containing only the Group1 columns for the specified omic modality, with column names renamed to match the original treatment name.
    result <- omicResults(simulation, omic_name) %>%
              select(!starts_with("Group2"))
    colnames(result) <- gsub("Group1.Time1", treatment, colnames(result))
    return(result)
}

simulation_one_donor <- function(rna_df, atac_df, IDtogene,
                                  outdir_simu_RNA, outdir_simu_ATAC,
                                  n_reps = 10) {
    #' Function to run the MOSim simulation for one donor by iterating over all
    #' treatments (columns) of the donor's data. 
    #' For each treatment it prepares the seed data by extracting the corresponding column from the RNA-seq and ATAC-seq data frames, 
    #' runs the MOSim simulation, post-processes the simulated data to select only the Group1 columns, 
    #' and saves the simulated data for each treatment in separate TSV files.
    #' @param rna_df Data frame of RNA-seq counts (genes x samples) for the specified donor.
    #' @param atac_df Data frame of ATAC-seq counts (peaks x samples) for the specified donor.
    #' @param IDtogene Data frame generated from the annotation file mapping peak IDs to gene IDs for the ATAC-seq data 
    #' @param outdir_simu_RNA Output directory for saving the simulated RNA-seq data for each treatment.
    #' @param outdir_simu_ATAC Output directory for saving the simulated ATAC-seq data for each treatment.
    #' @param n_reps Number of replicates to simulate for each treatment (default is 10).
    #' @return Invisible NULL. Writes simulated data files to the specified output directories for each treatment.

    for (current_ind_treat in colnames(atac_df)) {

        if (!current_ind_treat %in% colnames(rna_df)) {
            warning(paste("No:", current_ind_treat, "treatment found in the RNAseq data"))
            next
        }

        # Prepare seed data for both modalities
        rna_seed  <- prepare_seed_data(rna_df,  current_ind_treat)
        atac_seed <- prepare_seed_data(atac_df, current_ind_treat)
        cat("Library size RNA:", rna_seed$depth, "M | ATAC:", atac_seed$depth,
            "M for", current_ind_treat, "\n")

        # Run the MOSim simulation
        simulation <- mosim(
            omics        = c("RNA-seq", "DNase-seq"),
            times        = 1,
            numberGroups = 2,
            numberReps   = n_reps,
            diffGenes    = 0.0001,
            minMaxFC     = c(1.0, 1.0),
            omicsOptions = list(
                "RNA-seq"   = list(data = rna_seed$data,  depth = rna_seed$depth),
                "DNase-seq" = list(data = atac_seed$data, depth = atac_seed$depth,
                                   idToGene = IDtogene)
            )
        )

        # Post-process simulation results for both modalities
        dataRNAseq  <- postprocess_simulated_data(simulation, "RNA-seq",   current_ind_treat)
        dataATACseq <- postprocess_simulated_data(simulation, "DNase-seq", current_ind_treat)

        # Save the simulated data for this treatment
        write.table(dataRNAseq,
                    file.path(outdir_simu_RNA, paste0(current_ind_treat, "RNA.tsv")),
                    sep = "\t", quote = FALSE, col.names = NA)
        write.table(dataATACseq,
                    file.path(outdir_simu_ATAC, paste0(current_ind_treat, "ATAC.tsv")),
                    sep = "\t", quote = FALSE, col.names = NA)

        cat("Simulation of", current_ind_treat, "complete\n")
    }
    return(invisible(NULL))
}


merge_simulated_files <- function(sim_dir, individual_id, modality, output_dir) {
    #' Merge per-treatment simulated TSV files into a single file for one donor and modality.
    #' Verifies that all per-treatment files share the same feature IDs (genes or peaks)
    #' before concatenating them column-wise.
    #' @param sim_dir Directory containing the per-treatment simulated TSV files
    #'   (e.g., outdir_simu_RNA or outdir_simu_ATAC).
    #' @param individual_id Donor ID prefix to match in filenames (e.g., "HC19").
    #' @param modality Modality suffix in filenames: "RNA" or "ATAC".
    #' @param output_dir Directory where the merged TSV file will be saved.
    #' @return Invisible NULL. Writes a single merged TSV file to output_dir.
    #' @examples
    #' merge_simulated_files(outdir_simu_RNA, "HC19", "RNA", outdir_merged)

    # List per-treatment files for this donor and modality
    files <- list.files(sim_dir,
                        pattern = paste0("^", individual_id, ".*", modality, "\\.tsv$"),
                        full.names = TRUE)
    cat(modality, "files found:\n")
    cat(basename(files), sep = "\n")

    if (length(files) == 0) {
        stop("No ", modality, " files found in ", sim_dir, " for donor ", individual_id)
    }

    # Verify that all files share the same feature IDs (genes for RNA, peaks for ATAC)
    rownames_list <- lapply(files, function(f) {
        sort(rownames(read.table(f, header = TRUE, row.names = 1, sep = "\t")))
    })
    if (length(unique(rownames_list)) != 1) {
        stop("ERROR: ", modality, " files do not have the same feature IDs across treatments!")
    }
    
    # Concatenate column-wise after sorting rows for consistency
    merged <- do.call(cbind, lapply(files, function(f) {
        df <- read.table(f, header = TRUE, row.names = 1, sep = "\t")
        df[order(rownames(df)), , drop = FALSE]
    }))
    cat(modality, "merged:", nrow(merged), "x", ncol(merged), "\n")

    # Save the merged file
    output_path <- file.path(output_dir,
                             paste0(individual_id, "_", modality, "_merged.tsv"))
    write.table(merged, output_path,
                sep = "\t", quote = FALSE, col.names = NA)

    return(invisible(NULL))
}