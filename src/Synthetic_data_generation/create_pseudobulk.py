import argparse
import anndata
import numpy as np
import pandas as pd
import os

parser = argparse.ArgumentParser()
parser.add_argument("input_path", type=str)
parser.add_argument("output_path", type=str)
args = parser.parse_args()

os.makedirs(args.output_path, exist_ok=True)

path = args.input_path
filename = os.path.splitext(os.path.basename(path))[0]  # example "High_rep01_seed101_scDesign3"
parts = filename.split("_")
adnc = parts[0]   # "High"
rep  = parts[1]   # "rep01"
column_to_donor = f"{adnc}_{rep}"  # "High_rep01"
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
output_path = f"{args.output_path}/{adnc}_{rep}_pseudobulk.txt"

pseudobulk_df.to_csv(output_path, sep="\t")