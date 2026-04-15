"""Generate a pseudobulk count file from the simulated single-cell h5ad file for a given donor and replicate.

This script reads one simulated single-cell h5ad file, sums the counts across all cells for each gene to create a pseudobulk profile, 
and saves the result as a tab-separated file named after the donor and replicate.

Usage:
    python create_pseudobulk.py --input_path <input_path> --output_path <output_path>

Arguments:
    --input_path: Path to the input h5ad file. The filename must follow the
                  format <ADNC>_<rep>_<seed>_<tool>.h5ad (e.g. High_rep01_seed101_scDesign3.h5ad).
    --output_path: Directory where the pseudobulk .txt file will be saved.

Outputs:
    - <output_path>/<ADNC>_<rep>_pseudobulk.txt: Tab-separated pseudobulk count matrix
      with genes as rows and the sample label (<ADNC>_<rep>) as column header.
"""


import argparse
import anndata
import numpy as np
import pandas as pd
import os


parser = argparse.ArgumentParser(
    description="Generate a pseudobulk count file from a simulated donor h5ad file."
)
parser.add_argument("--input_path", type=str,
    help="Path to the input h5ad file. Filename must follow <ADNC>_<rep>_<seed>_<tool>.h5ad format.")
parser.add_argument("--output_path", type=str,
    help="Directory where the pseudobulk .txt file will be saved.")
args = parser.parse_args()

# Create the output directory if it doesn't exist
os.makedirs(args.output_path, exist_ok=True)

# Extract the donor name and replicate number from the input file name
path = args.input_path
filename = os.path.splitext(os.path.basename(path))[0]  # example "High_rep01_seed101_scDesign3"
parts = filename.split("_")
adnc = parts[0]   # "High"
parts = filename.split("_")
adnc = parts[0]
print(f"Processing: {filename}, parts: {parts}")
rep = next(p for p in parts if p.startswith("rep"))  # "rep01"
column_to_donor = f"{adnc}_{rep}"  # "High_rep01"

# Load the single-cell data from the input file
adata = anndata.read_h5ad(path)

# calculate pseudobulk counts by summing the counts across all cells for each gene
pseudobulk = np.array(adata.X.sum(axis=0)).flatten()

#create a DataFrame with the pseudobulk counts
pseudobulk_df = pd.DataFrame(
    pseudobulk,
    index=adata.var_names,
    columns=[column_to_donor],
    dtype=np.int64
)

# Save the pseudobulk results to a text file
output_path = os.path.join(args.output_path, f"{adnc}_{rep}_pseudobulk.txt")
pseudobulk_df.to_csv(output_path, sep="\t")  