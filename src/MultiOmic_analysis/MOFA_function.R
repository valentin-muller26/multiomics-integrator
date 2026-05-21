train_mofa_model  <- function(omics_df, out_dir, nb_factor = 10, seed = 42, nb_iter = 1000){
  #' Function for training the MOFA model
  #'@param omics_df : long format dataframe containing all the omics data
  #'@param out_dir : output directory for the model fitted in format hdf5
  #'@param nb_factor : number of factor that the mofa models will used to explained the variance in the data default [default=10]
  #'@param seed : seed used for reproducibility of the fitting model [default = 42]
  #'@param nb_iter : number of iteration in the fitting step [default = 1000]
  #'@return MOFAobject.trained : return the trained model 
  
  #Create the mofa object
  MOFAobject <- create_mofa(omics_df)
  
  #Use the default option for the model execpt scale_views since the two omics used different normalisation method
  data_opts <- get_default_data_options(MOFAobject)
  data_opts$scale_views  <- TRUE     # VST vs log-CPM
  data_opts$scale_groups <- FALSE
  model_opts <- get_default_model_options(MOFAobject)
  model_opts$num_factors <- nb_factor
  
  # Use default parameter for the model training and set the convergence mode to slow and the number of iteration to the one set by the user
  train_opts <- get_default_training_options(MOFAobject)
  train_opts$seed             <- seed
  train_opts$convergence_mode <- "slow"
  train_opts$maxiter          <- nb_iter
  
  #Object for the parameter of the training model 
  MOFAobject <- prepare_mofa(
    object           = MOFAobject,
    data_options     = data_opts,
    model_options    = model_opts,
    training_options = train_opts
  )
  
  # create the output directory if not present
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  outfile <- file.path(out_dir, "model_donor_corrected.hdf5")
  #Run the training of the mode 
  MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk = TRUE)
  return(MOFAobject.trained)
}

add_metadata <- function(MOFAobject.trained){
  #' Add sample metadata to a trained MOFA model
  #'
  #' Extracts donor identifiers (HC followed by digits) and cytokine condition
  #' from sample names, then attaches them as metadata to the MOFA object.
  #' The `group` column is set to "single_group" because the model is trained
  #' on a single group (no group-wise factorization).

  #'@param MOFAobject.trained : A trained MOFA2 model object whose sample names follow the pattern HC<number of patient>_<treatment>
  #'@return MOFAobject.trained : trained model with the metadata
  
  #Retrieve the sample name from the trained mdel 
  sample_names <- unlist(samples_names(MOFAobject.trained)) 

  # Build the metadata table by parsing donor and condition from sample names to guarante the same order
  samples_meta <- data.frame(
    sample    = sample_names,
    group     = "single_group",
    donor     = str_extract(sample_names, "HC[0-9]+"),
    condition = str_remove(sample_names, "HC[0-9]+")
  )
  samples_metadata(MOFAobject.trained) <- samples_meta
  return(MOFAobject.trained)
}


extract_pathway_genes <- function(enrichment, factor_name, 
                                  n_top_pathways = 25, 
                                  n_top_genes_per_pathway = 50) {

  #' Extract top pathways and their top genes from a MOFA2 GSEA enrichment result
  #'@description
  #' For a given MOFA2 factor, selects the most significant pathways ranked by
  #' adjusted p-value and, within each pathway, the genes contributing most
  #' strongly to the factor ranked by absolute feature statistic. Each gene is
  #' annotated with its gene name symbol and description via org.Hs.eg.db when an
  #' Ensembl ID match is available, otherwise the original ID is kept.
  #'
  #'@param enrichment : MOFA2 enrichment result returned by run_enrichment()
  #'@param factor_name : name of the MOFA2 factor to summarise (e.g. "Factor1")
  #'@param n_top_pathways : maximum number of pathways to keep, ordered by ascending adjusted p-value [default = 25]
  #'@param n_top_genes_per_pathway : maximum number of genes to keep per pathway, ordered by descending absolute feature statistic [default = 50]
  #'@return pathway_table : data.frame with one row per pathway containing the pathway name, raw and adjusted p-values, set statistic, number of genes, and the top genes as Ensembl IDs, gene name and description
  #---------Extracting the component for the factor-------------------------------
  
  
  #Extract matrix indicating if the gene are part of the pathway (0 absent 1 present)
  pathway_gene_membership <- enrichment$feature.sets
  
  # extract p-value adjust and value for the enrichisment for each pathway for the factor
  pathway_padj <- enrichment$pval.adj[, factor_name]
  pathway_test_statistic <- enrichment$set.statistics[, factor_name]
  gene_factor_statistic <- enrichment$feature.statistics[, factor_name]
  
  # ---- Select the n pathway most significative ----------------------------------------------------
  
  pathway_order_by_significance <- order(pathway_padj, decreasing = FALSE)
  
  # Keep n pathway and if less take all the pathway
  n_pathways_to_keep <- min(n_top_pathways, length(pathway_order_by_significance))
  
  # retrieve name of the pathway by order of p_adj value and select the top n pathway
  top_pathway_names <- rownames(pathway_gene_membership)[pathway_order_by_significance]
  top_pathway_names <- top_pathway_names[seq_len(n_pathways_to_keep)]
  
  # ---- Construct the taabel ----
  
  # list to accumulate the data (one per line/pathway)
  pathway_rows <- list()
  
  for (pathway_name in top_pathway_names) {
    
    # Retrieve the gene present in the pathway
    is_gene_in_pathway   <- pathway_gene_membership[pathway_name, ] == 1
    genes_in_pathway     <- colnames(pathway_gene_membership)[is_gene_in_pathway]
    n_genes_in_pathway   <- length(genes_in_pathway)
    
    if (n_genes_in_pathway > 0) {
      
      # Get the statistic of the gene in pathway
      gene_stats_in_pathway <- gene_factor_statistic[genes_in_pathway]
      
      # Suppress the Na 
      gene_stats_in_pathway <- gene_stats_in_pathway[!is.na(gene_stats_in_pathway)]
      
      # Sort the gene by statistic to have the gene that contribute the most first
      gene_order_by_strength <- order(abs(gene_stats_in_pathway), decreasing = TRUE)
      
      # Keep the n top gene
      n_genes_to_keep <- min(n_top_genes_per_pathway, length(gene_stats_in_pathway))
      top_gene_names  <- names(gene_stats_in_pathway)[gene_order_by_strength]
      top_gene_names  <- top_gene_names[seq_len(n_genes_to_keep)]
      
      annot <- suppressMessages(AnnotationDbi::select(
        org.Hs.eg.db,
        keys    = top_gene_names,
        columns = c("SYMBOL", "GENENAME"),
        keytype = "ENSEMBL"
      ))
      annot <- annot[!duplicated(annot$ENSEMBL), ]
      rownames(annot) <- annot$ENSEMBL
      
      top_gene_symbols   <- annot[top_gene_names, "SYMBOL"]
      top_gene_fullnames <- annot[top_gene_names, "GENENAME"]
      top_gene_symbols[is.na(top_gene_symbols)]     <- top_gene_names[is.na(top_gene_symbols)]
      top_gene_fullnames[is.na(top_gene_fullnames)] <- ""
      
      top_genes_concatenated        <- paste(top_gene_names,     collapse = ", ")
      top_genes_symbol_concatenated <- paste(top_gene_symbols,   collapse = ", ")
      top_genes_name_concatenated   <- paste(top_gene_fullnames, collapse = "; ")
      
    } else {
      # Case pathway with no gene
      top_genes_concatenated <- ""
      top_genes_symbol_concatenated <- ""
      top_genes_name_concatenated   <- ""
    }
    
    # Construct the dataframe
    pathway_rows[[pathway_name]] <- data.frame(
      pathway        = pathway_name,
      pval           = enrichment$pval[pathway_name, factor_name],
      pval_adj       = pathway_padj[pathway_name],
      set_statistic  = pathway_test_statistic[pathway_name],
      n_genes_total  = n_genes_in_pathway,
      top_genes      = top_genes_concatenated,
      top_genes_symbol = top_genes_symbol_concatenated,
      top_genes_name   = top_genes_name_concatenated,
      stringsAsFactors = FALSE
    )
  }
  
  # bind each row in one dataframe table
  pathway_table <- do.call(rbind, pathway_rows)
  return(pathway_table)
}