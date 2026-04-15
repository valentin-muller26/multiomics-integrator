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

os.makedirs(args.output_path, exist_ok=True)

def read_and_merge(directory):
    files = sorted(
        os.path.join(directory, f)
        for f in os.listdir(directory)
        if f.endswith(".txt")
    )
    dfs = [pd.read_csv(f, sep="\t", index_col=0) for f in files]

    ref_index = dfs[0].index
    for f, df in zip(files[1:], dfs[1:]):
        if not ref_index.equals(df.index):
            only_ref = len(ref_index.difference(df.index))
            only_other = len(df.index.difference(ref_index))
            print(f"WARNING: index mismatch with {os.path.basename(f)} "
                  f"(+{only_ref} only in ref, +{only_other} only in this file) → outer join + fillna(0)")

    merged = pd.concat(dfs, axis=1, join="outer").fillna(0).astype(int)
    merged = merged.sort_index(axis=0)  # consistent feature order
    return merged

merged_rna  = read_and_merge(args.input_path_rna)
merged_atac = read_and_merge(args.input_path_atac)

merged_rna.to_csv( os.path.join(args.output_path, "real_RNA_merged_pseudobulk.txt"),  sep="\t")
merged_atac.to_csv(os.path.join(args.output_path, "real_ATAC_merged_pseudobulk.txt"), sep="\t")

print(f"RNA : {merged_rna.shape[1]} samples x {merged_rna.shape[0]} genes")
print(f"ATAC: {merged_atac.shape[1]} samples x {merged_atac.shape[0]} peaks")