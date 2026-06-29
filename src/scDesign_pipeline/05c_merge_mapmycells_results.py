"""
Lists all CSV files in the specified input directory that match the pattern "*_scDesign3_percentage_per_celltype.csv",
 merges them into a single DataFrame, and saves the merged results to a new CSV containing the percentage of correctly classified cells per cell type for each replicate.
Additionally, this script calculates summary statistics (min, max, mean, variance) across replicates for each cell type and saves these statistics 
to a separate CSV file in the output directory.


Usage:
    python 05c_merge_mapmycells_results.py --input <input_directory> --output <output_directory>
"""
import argparse
import pandas as pd
import os
import glob

parser = argparse.ArgumentParser()
parser.add_argument("--input", help="Input directory containing MapMyCells results")
parser.add_argument("--output", help="Output directory for merged results")
args = parser.parse_args()


# Create output directory if it doesn't exist
os.makedirs(args.output, exist_ok=True)
# Read all CSV files in the input directory that match the pattern and merge them into a single DataFrame
files = sorted(glob.glob(os.path.join(args.input, "*_scDesign3_percentage_per_celltype.csv")))

reference_order = None
columns = {}

for file in files:
    # Read the CSV file and sort by cell_type to ensure consistent order
    df = pd.read_csv(file).sort_values("cell_type").reset_index(drop=True)
    # Extract replicate name from the filename
    replicate = os.path.basename(file).split("_scDesign3_percentage_per_celltype.csv")[0]
    ADNC_group = replicate.split("_")[0] 
    # Verify cell types are in the same order across all files
    if reference_order is None:
        reference_order = df["cell_type"].tolist()

    columns[replicate] = df["percentage_correct"].values

df_merged = pd.DataFrame({"cell_type": reference_order})
for replicate, values in columns.items():
    df_merged[replicate] = values

df_merged.to_csv(os.path.join(args.output, f"{ADNC_group}_merged_scDesign3_percentage_per_celltype.csv"), index=False)

# Calculate summary statistics (min, max, mean, variance) across replicates for each cell type
replicate_cols = [c for c in df_merged.columns if c != "cell_type"]

df_stat = pd.DataFrame({"cell_type": reference_order})
df_stat["min"] = df_merged[replicate_cols].min(axis=1)
df_stat["max"] = df_merged[replicate_cols].max(axis=1)
df_stat["mean"] = df_merged[replicate_cols].mean(axis=1)
df_stat["variance"] = df_merged[replicate_cols].var(axis=1)


df_stat.to_csv(os.path.join(args.output, f"{ADNC_group}_summary_statistics_mapmycells.csv"), index=False)