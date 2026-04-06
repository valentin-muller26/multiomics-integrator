"""Visualize cell type distribution across donors as a clustered heatmap.

This script reads the individual donor RNA h5ad files, extract the cell type distribution (Subclass) and the ADNC status for each donor, 
and then creates a clustered heatmap to visualize the cell type distribution across donors, with ADNC status as column colors, 
the color of each cell representing the percentage of cell types for the specific donor and raw cellcounts indicated inside each cell of the heatmap.

Usage:
    python 02a_donor_cell_distribution.py <input_folder> <output_folder>

Arguments:
    input_folder: Directory containing per-donor *_RNA_multiome_subset.h5ad files.
    output_folder: Directory where the output heatmap PNG will be saved.

Outputs:
    - <output_folder>/celltype_distribution.png: Clustered heatmap of normalized
      cell type proportions (%) with raw counts as values inside each cell of the heatmap.
"""

import anndata
import pandas as pd
from matplotlib.lines import Line2D
import seaborn as sns
import os
import argparse



def analyze_donor(file_path):
    """Extract cell type distribution and metadata from a donor h5ad file.

    Args:
        file_path: Path to the donor RNA h5ad file.

    Returns:
        A tuple of:
            - subclass_counts (pd.Series): Cell counts per Subclass /celltype label.
            - donor_id (str): Unique donor identifier.
            - adnc (str): Overall AD neuropathological Change category.
    """
    a = anndata.read_h5ad(file_path)
    donor_id = a.obs['Donor ID'].unique()[0]
    ADNC = a.obs["Overall AD neuropathological Change"].unique()[0]
    return a.obs["Subclass"].value_counts(), donor_id,ADNC

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot clustered heatmap of cell type distribution per donor."
    )
    parser.add_argument("input_folder", type=str, help="Folder with per-donor RNA h5ad files.")
    parser.add_argument("output_folder", type=str, help="Folder to save the output heatmap.")
    args = parser.parse_args()

    # List all the donor subset RNA files in the input folder
    file_paths = [os.path.join(args.input_folder, f) for f in os.listdir(args.input_folder) if f.endswith("RNA_multiome_subset.h5ad")]


    #Temporary dictionaries to store the cell type distribution and ADNC status for each patient
    donor_ADNC_status = dict()
    donor_cell_types = dict()

    # Analyze each donor file to extract cell type distribution and ADNC status and donor ID
    for donor_path in file_paths : 
        subclass,donor_ID,ADNC= analyze_donor(donor_path)
        donor_cell_types[donor_ID] = subclass
        donor_ADNC_status[donor_ID] = ADNC

    # Convert the dictionaries to DataFrames for easier manipulation and plotting
    donor_cell_types = pd.DataFrame(donor_cell_types)
    donor_ADNC_status  = pd.Series(donor_ADNC_status)


    # Define a color palette for the ADNC categories
    ADNC_COLORS = {
        "Not AD":       "#2196F3",  # bleu
        "Low":          "#4CAF50",  # vert
        "Intermediate": "#FF9800",  # orange
        "High":         "#F44336",  # rouge
    }


    #Mapping the ADNC status to colors for the heatmap annotation
    col_colors = donor_ADNC_status.map(ADNC_COLORS)

    # Normalize the cell type counts to percentages for better visualization in the heatmap and fill NaN values with 0 and force type to int
    donor_cell_types.fillna(0, inplace=True)
    donor_cell_types = donor_cell_types.astype(int)
    donor_cell_types_norm = donor_cell_types.div(donor_cell_types.sum(axis=0), axis=1) * 100

    # Create a clustered heatmap to visualize the cell type distribution across patients, with ADNC status as column colors and y axis labels as cell types and x axis labels as donor  IDs with total number of cells in each patient
    heatmap = sns.clustermap(donor_cell_types_norm, 
                   col_colors=col_colors, 
                   annot=donor_cell_types, 
                   fmt='d',
                   cmap="coolwarm",
                   figsize=(20, 12)) 

    #Reorder the columns based on the clustering results to display the cell type distribution in the same order as the heatmap
    reordered_cols = donor_cell_types_norm.columns[heatmap.dendrogram_col.reordered_ind]
    # Version not normalized
    # heatmap = sns.clustermap(donor_cell_types, 
    #                col_colors=col_colors, 
    #                annot=donor_cell_types, 
    #                fmt='d',
    #                cmap="coolwarm",
    #                figsize=(20, 12)) 
    # reordered_cols = donor_cell_types.columns[heatmap.dendrogram_col.reordered_ind]


    # Update the x-axis labels to include the number of cells for each patient and customize the label
    new_labels = [f"{col}\n(n={donor_cell_types[col].sum()})" for col in reordered_cols]
    heatmap.ax_heatmap.set_xticklabels(new_labels, rotation=45, ha='right', fontsize=7)


    # Create custom legend handles for the ADNC categories
    legend_handles = [
        Line2D([0], [0], marker='s', color='w',
            markerfacecolor=color, markersize=11, label=label)
        for label, color in ADNC_COLORS.items()
    ]

    # Add the legend to the heatmap
    heatmap.ax_col_dendrogram.legend(
    handles=legend_handles,
    title="ADNC",
    loc="center",
    frameon=False,
    )
    # Save the heatmap to the output folder created if it doesn't exist
    os.makedirs(args.output_folder, exist_ok=True)

    heatmap.savefig(os.path.join(args.output_folder, "celltype_distribution.png"), dpi=150, bbox_inches="tight")