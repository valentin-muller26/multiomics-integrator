"""Merge pseudobulk RNA and ATAC files into a single DataFrame and generate metadata.

This script reads individual pseudobulk RNA and ATAC .txt files from RNA and ATAC input directories, 
merges them into a single files for the corresponding modality and reorders the ATAC columns to match the RNA sample order.
It also generates a metadata file with the ADNC condition and numeric label for each sample, which can be used for downstream analyses and deep learning models.

Usage:
    python merge_pseudobulk.py <input_path_rna> <input_path_atac> <output_path>

Arguments:
    input_path_rna: Directory containing per-donor-replicate pseudobulk RNA .txt files.
    input_path_atac: Directory containing per-donor-replicate pseudobulk ATAC .txt files.
    output_path: Directory where the merged files and metadata will be saved.

Outputs:
    - <output_path>/RNA_merged_pseudobulk.txt: Merged RNA pseudobulk matrix.
    - <output_path>/ATAC_merged_pseudobulk.txt: Merged ATAC pseudobulk matrix.
    - <output_path>/metadata.csv: Sample metadata with ADNC condition and numeric label.
"""


import argparse
import pandas as pd
import os


ADNC_LABELS = {
    "Not AD": 0,
    "Low": 1,
    "Intermediate": 2,
    "High": 3,
}


def load_and_merge(input_path: str) -> pd.DataFrame:
    """Read all .txt pseudobulk files in input_path and merge them into a single DataFrame.

    Each file is expected to be a tab-separated matrix with features as rows and one sample as column.
    Rows are sorted before merging to ensure consistent index alignment across files.

    Args:
        input_path: Directory containing per-donor-replicate pseudobulk .txt files.

    Returns:
        Merged DataFrame with features as rows and samples as columns, sorted by column name.
    """
    files = [
        os.path.join(input_path, f)
        for f in os.listdir(input_path)
        if f.endswith(".txt")
    ]
    dfs = []
    for file in files:
        df = pd.read_csv(file, sep="\t", index_col=0)
        df = df.sort_index()
        dfs.append(df)

    merged = pd.concat(dfs, axis=1, join="outer").fillna(0).astype(int)
    return merged.sort_index(axis=1)


parser = argparse.ArgumentParser(
    description="Merge pseudobulk RNA and ATAC files and generate metadata."
)
parser.add_argument("input_path_rna", type=str, help="Directory containing per-donor-replicate pseudobulk RNA .txt files.")
parser.add_argument("input_path_atac", type=str, help="Directory containing per-donor-replicate pseudobulk ATAC .txt files.")
parser.add_argument("output_path", type=str, help="Directory where the merged files and metadata will be saved.")
args = parser.parse_args()

os.makedirs(args.output_path, exist_ok=True)

# --- 1. Merge RNA and ATAC pseudobulk files ---
merged_rna = load_and_merge(args.input_path_rna)
merged_atac = load_and_merge(args.input_path_atac)

# --- 2. Check and align samples across modalities ---
rna_samples = set(merged_rna.columns)
atac_samples = set(merged_atac.columns)

missing_in_atac = rna_samples - atac_samples
missing_in_rna = atac_samples - rna_samples

if missing_in_atac:
    print(f"WARNING: {len(missing_in_atac)} RNA samples missing in ATAC: {missing_in_atac}")
if missing_in_rna:
    print(f"WARNING: {len(missing_in_rna)} ATAC samples missing in RNA: {missing_in_rna}")

# Keep only common samples and reorder ATAC columns to match RNA sample order
common_samples = sorted(rna_samples & atac_samples)
merged_rna = merged_rna[common_samples]
merged_atac = merged_atac[common_samples]

# Save merged files
merged_rna.to_csv(os.path.join(args.output_path, "RNA_merged_pseudobulk.txt"), sep="\t")
merged_atac.to_csv(os.path.join(args.output_path, "ATAC_merged_pseudobulk.txt"), sep="\t")

# --- 3. Metadata ---
metadata = pd.DataFrame({
    "donor_rep": common_samples,
    "ADNC": [col.split("_")[0] for col in common_samples],
})
metadata["label"] = metadata["ADNC"].map(ADNC_LABELS)
metadata.to_csv(os.path.join(args.output_path, "metadata.csv"), index=False)