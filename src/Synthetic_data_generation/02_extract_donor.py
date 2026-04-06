"""Extract and save paired Multiome RNA and ATAC subsets per donor from SEA-AD dataset.

This script processes a pair of RNA and ATAC h5ad files from the SEA-AD dataset containing multi donors and multiple modalities (10x Multiome and snRNA-seq or ATACseq data)
This script identifies the Multiome donors subsets in RNAseq h5ad and then extracts the paired barcodes in the ATAC h5ad, ensuring that only cells profiled in both modalities are retained.
For each donor, it saves the paired RNA and ATAC subsets as separate h5ad files. ¨
Additionally, it compiles a summary CSV file containing metadata for all processed donors and a flat list of donor IDs for downstream SLURM array jobs.

Usage:
    python 02_extract_donor.py <input_rna_path> <input_atac_path> <output_path>

Args:
    input_rna_path: Path to the RNA h5ad file.
    input_atac_path: Path to the ATAC h5ad file.
    output_path: Directory where per-donor h5ad files and summary files(metadata and donor list) will be saved.

Outputs:
    - <output_path>/<donor>_RNA_multiome_subset.h5ad: Per-donor RNA subset.
    - <output_path>/<donor>_ATAC_multiome_subset.h5ad: Per-donor ATAC subset.
    - <output_path>/donors_summary.csv: Summary metadata for all processed donors.
    - <output_path>/donor_ids.txt: Flat list of donor IDs for SLURM array jobs.
"""


import argparse
import anndata as ad
import pandas as pd
import gc
from pathlib import Path

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_rna_path", type=str)
    parser.add_argument("input_atac_path", type=str)
    parser.add_argument("output_path", type=str)
    args = parser.parse_args()

    # Create the output directory if it doesn't exist
    output_dir = Path(args.output_path)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load the RNA and ATAC data in backed mode to handle large files efficiently
    print("Loading in backed mode...")
    adata_rna = ad.read_h5ad(args.input_rna_path, backed='r')
    adata_atac = ad.read_h5ad(args.input_atac_path, backed='r')

    # Summarize of the loaded data
    print(f"RNA cells: {adata_rna.n_obs}")
    print(f"ATAC cells: {adata_atac.n_obs}")
    print(f"Available Methods: {adata_rna.obs['method'].unique().tolist()}")
    print(f"Available Donors: {adata_rna.obs['Donor ID'].unique().tolist()}")

    # Identify Multiome cells in the RNA data
    multiome_mask = adata_rna.obs['method'] == '10xMulti'
    multiome_barcodes = adata_rna.obs_names[multiome_mask]
    print(f"Multiome cells in RNA: {len(multiome_barcodes)}")

    # Intersections of barcodes between RNA and ATAC to ensure we only keep paired cells
    paired_barcodes = multiome_barcodes.intersection(adata_atac.obs_names)
    print(f"Multiome cells paired in ATAC: {len(paired_barcodes)}")


    if len(paired_barcodes) == 0:
        print("ERROR: No Multiome barcodes found in ATAC file.")
        print(f"RNA Multiome obs_names examples: {multiome_barcodes[:5].tolist()}")
        print(f"ATAC obs_names examples: {adata_atac.obs_names[:5].tolist()}")
        exit(1)

    # Mask final : Multiome + present in both RNA and ATAC
    paired_mask = multiome_mask & adata_rna.obs_names.isin(paired_barcodes)

    multiome_donors = adata_rna.obs.loc[paired_mask, 'Donor ID'].unique()
    print(f"Multiome donors with paired data: {len(multiome_donors)}")

    # Collecting information for each donor and saving subsets
    optional_cols = [
        'Age at Death', 'Sex', 'method', 'Overall AD neuropathological Change',
        'Braak', 'CERAD score', 'Cognitive Status', 'Neurotypical reference',
        'Class', 'Subclass', 'Supertype', 'Brain Region', 'APOE Genotype',
        'Continuous Pseudo-progression Score', 'Severely Affected Donor'
    ]
    list_donors_information = []
    # Iterate over each donor and extract the paired RNA and ATAC subsets, and metadata information and save them
    for donor in multiome_donors:
        combined_mask = paired_mask & (adata_rna.obs['Donor ID'] == donor)
        n_cells = combined_mask.sum()
        print(f"\nDonor: {donor}, Cells: {n_cells}")

        if n_cells == 0:
            print(f"  Skipping {donor}: no paired cells.")
            continue

        # Extraction RNA
        donor_rna_subset = adata_rna[combined_mask].to_memory()

        # Extraction ATAC via obs_names
        donor_barcodes = donor_rna_subset.obs_names
        donor_atac_subset = adata_atac[donor_barcodes].to_memory()

        # verify that the number of cells match in both subsets
        assert donor_rna_subset.n_obs == donor_atac_subset.n_obs, (
            f"Mismatch for {donor}: RNA={donor_rna_subset.n_obs}, ATAC={donor_atac_subset.n_obs}"
        )
        print(f"  Paired cells verified: {donor_rna_subset.n_obs}")

        # Collecting main donor information
        donor_info = {
            'Donor ID': donor,
            'cell_count': int(combined_mask.sum()),
            'n_genes': donor_rna_subset.n_vars,
            'n_peaks': donor_atac_subset.n_vars,
        }
        # Adding optional metadata if available
        for col in optional_cols:
            if col in donor_rna_subset.obs.columns:
                vals = donor_rna_subset.obs[col].unique().tolist()
                donor_info[col] = "; ".join(str(v) for v in vals)
            else:
                donor_info[col] = "N/A"

        list_donors_information.append(donor_info)

        # Save the donor RNA subset
        out_rna_path = output_dir / f"{donor}_RNA_multiome_subset.h5ad"
        donor_rna_subset.write_h5ad(out_rna_path)
        print(f"Saved RNA: {out_rna_path}")
        # Save the donor ATAC subset
        out_atac_path = output_dir / f"{donor}_ATAC_multiome_subset.h5ad"
        donor_atac_subset.write_h5ad(out_atac_path)
        print(f"Saved ATAC: {out_atac_path}")

        # Force garbage collection to free memory
        del donor_rna_subset, donor_atac_subset
        gc.collect()

    # Save the summary of donors information
    df_summary = pd.DataFrame(list_donors_information)
    data_summary_path = output_dir / "donors_summary.csv"
    df_summary.to_csv(data_summary_path, index=False)
    print(f"Summary saved: {data_summary_path}")

    # Save the list of donor IDs for array jobs in the next steps
    df_summary['Donor ID'].to_csv(output_dir / "donor_ids.txt", index=False, header=False)
    print(f"List of donor IDs saved: {output_dir / 'donor_ids.txt'}")