"""Generate a pseudobulk count file from a real single-cell h5ad file.

Reads one real h5ad file (used as simulation seed), sums counts across all cells
per gene and divides by the number of donors to normalize, and saves the result as a tab-separated pseudobulk profile.

Usage:
    python create_pseudobulk_real.py --input_path <input_path> --output_path <output_path> --modality <RNA|ATAC>

Arguments:
    --input_path:  Path to the input h5ad file. Filename must follow real_<ADNC>.h5ad
                   (e.g. real_High.h5ad, real_Not_AD.h5ad).
    --output_path: Directory where the pseudobulk .txt file will be saved.
    --modality:    Modality: RNA (uses layers['UMIs']) or ATAC (uses X).

Outputs:
    - <output_path>/real_<ADNC>_pseudobulk.txt
"""

import argparse
import anndata
import numpy as np
import pandas as pd
import os

parser = argparse.ArgumentParser(description="Generate a pseudobulk count file from a real donor h5ad file.")
parser.add_argument("--input_path", type=str,
    help="Path to the input h5ad file. Filename must follow real_<ADNC>.h5ad format.",required=True)
parser.add_argument("--output_path", type=str,
    help="Directory where the pseudobulk .txt file will be saved.",required=True)
parser.add_argument("--modality", type=str, choices=["RNA", "ATAC"], required=True,
    help="Modality: RNA (uses layers['UMIs']) or ATAC (uses X).")
args = parser.parse_args()


# Parse filename: real_High.h5ad -> label = "real_High"
filename = os.path.splitext(os.path.basename(args.input_path))[0]  # "real_High"
assert filename.startswith("real_"), f"Unexpected filename format: {filename}"
label = filename  # keep "real_<ADNC>" as column name
print(f"Processing: {filename}")

# Read the h5ad file
adata = anndata.read_h5ad(args.input_path)
# For RNA, use the "UMIs" layer which contains the raw counts. For ATAC, use the X matrix which contains the raw counts.
if args.modality == "RNA":
    mat = adata.layers["UMIs"]
else:
    mat = adata.X

# Number of unique donors in this ADNC group
n_donors = adata.obs["Donor.ID"].nunique()
print(f"Number of donors in {label}: {n_donors}")

# Sum counts across all cells per gene and divide by the number of donors to normalize
pseudobulk = np.array(mat.sum(axis=0)).flatten()
pseudobulk_normalized = np.round(pseudobulk / n_donors).astype(np.int64)

pseudobulk_df = pd.DataFrame(
    pseudobulk_normalized,
    index=adata.var_names,
    columns=[label],
    dtype=np.int64
)

# Save the pseudobulk as a tab-separated .txt file
# Create output directory if it doesn't exist
os.makedirs(args.output_path, exist_ok=True)
output_file = os.path.join(args.output_path, f"{filename}_pseudobulk.txt")
pseudobulk_df.to_csv(output_file, sep="\t")
print(f"Saved: {output_file}")