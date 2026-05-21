"""Merge real pseudobulk RNA and ATAC files into a single DataFrame each.

Usage:
    python merge_pseudobulk_real.py <input_path_rna> <input_path_atac> <output_path>

Outputs:
    - <output_path>/real_RNA_merged_pseudobulk.txt
    - <output_path>/real_ATAC_merged_pseudobulk.txt
"""

import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser()
parser.add_argument("input_path_rna", type=str)
parser.add_argument("input_path_atac", type=str)
parser.add_argument("output_path", type=str)
args = parser.parse_args()


def read_and_merge(directory):
    """Read all the pseudobulk .txt files in the given directory, merge them into a single DataFrame, and return it.
        use the first file as reference for the features name and order, and compare the others to it. 
        If there are mismatches in the features, print a warning and perform an outer join with fillna(0) to include all features.
    Args :
        directory (str): Path to the directory containing the pseudobulk .txt files.
    Usage:
        read_and_merge("/path/to/pseudobulk/files")
    Outputs:
        A merged DataFrame containing all samples from the .txt files, with features as rows ordered alphabetically to have a consistent order 
        with the merged pseudobulk simulated data and samples as columns 
    """
    # Read all .txt files in the directory into a list of DataFrames, ensuring they are sorted for consistent order
    files = sorted(
        os.path.join(directory, f)
        for f in os.listdir(directory)
        if f.endswith(".txt")
    )
    if not files:
        raise FileNotFoundError(f"No .txt files found in {directory}")
    dfs = [pd.read_csv(f, sep="\t", index_col=0) for f in files]


    # Check if all DataFrames have the same index (features). If not, print a warning and perform an outer join with fillna(0) to include all features.
    # use the first file as reference for the features name and order, and compare the others to it
    ref_index = dfs[0].index
    for file, df in zip(files[1:], dfs[1:]):
        if not ref_index.equals(df.index):
            only_ref = len(ref_index.difference(df.index))
            only_other = len(df.index.difference(ref_index))
            print(f"WARNING: index mismatch with {os.path.basename(file)} "
                  f"(+{only_ref} only in ref, +{only_other} only in this file) → outer join + fillna(0)")

    merged = pd.concat(dfs, axis=1, join="outer").fillna(0).astype(int) # convert to int after filling with 0 the absence of features in some files
    merged = merged.sort_index(axis=0)  # sort features alphabetically for consistency 
    return merged

# Read and merge RNA and ATAC pseudobulk files
merged_rna  = read_and_merge(args.input_path_rna)
merged_atac = read_and_merge(args.input_path_atac)

# Save merged matrices to output directory creating it if it doesn't exist
os.makedirs(args.output_path, exist_ok=True)
merged_rna.to_csv( os.path.join(args.output_path, "real_RNA_merged_pseudobulk.txt"),  sep="\t")
merged_atac.to_csv(os.path.join(args.output_path, "real_ATAC_merged_pseudobulk.txt"), sep="\t")

print(f"RNA : {merged_rna.shape[0]} genes x {merged_rna.shape[1]} samples")
print(f"ATAC: {merged_atac.shape[0]} peaks x {merged_atac.shape[1]} samples")