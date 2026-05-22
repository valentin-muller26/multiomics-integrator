# MultiOmic analysis
This repository consists of the  Code and configuration for three pipelines developed as part of a Master's thesis in Bioinformatics and Computational Biology (Universities of Fribourg & Berne) on the multi-omics integration of bulk RNA-seq and ATAC-seq data from healthy donors in the context of atopic dermatitis.

## Repository overview

The repository is divided in three independant pipeline : 
1.**Single-cell dataset simulation** : in this pipeline, the `SEA-AD MTG Multiome` dataset was used as a proof-of-concept dataset to evaluate the feasibility of simulating paired single-cell RNA-seq and ATAC-seq data with scDesign3 and consiss of the following step : downloading, preprocessing and quality control of the SEA-AD MTG Multiome dataset, simulation with `scDesign3`, and validation of the simulated counts.
2. **Simulation of the bulk RNA-seq and ATAC-seq data from healthy donors dataset** : Preprocessing of the bulk RNA-seq and ATAC-seq from healthy donors in the context of atopic dermatitis, simulation with `MOSim`, and validation with `countsimQC`. 
3. **Multi-omics analysis of the bulk RNA-seq and ATAC-seq data from healthy donors dataset** :  Preprocessing of the bulk RNA-seq and ATAC-seq from healthy donors, training and downstream analysis of  unsupervised `MOFA2` models, and training of the supervised model `DIABLO` (mixOmics).



Pipeline of single cell dataset simulation consisting of downloading the dataset of SEA-AD, preprocessing the dataset and quality control, simulation using scDesign3 and validation of the simulation.
3. Pipeline of the simulation of the dataset containing of omics (RNA-seq and ATAC-seq) of healthy donor for atopic dermatitis consisting of the step of preprocessing, simulation using MOSim, and validation using countSimQC.
4. Pipeline for the Muliomics analysis of the dataset containing of omics (RNA-seq and ATAC-seq) of healthy donor for atopic dermatitis consisting of the step of preprocessing the datasset, training the model of MOFA and downstream analysis of MOFA, and analysis using Diablo.

## Structure of the repository 

