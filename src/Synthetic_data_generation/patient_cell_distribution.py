import anndata
import pandas as pd
from matplotlib.lines import Line2D
import seaborn as sns
import os
import argparse



def analyze_patient(file_path):
    a = anndata.read_h5ad(file_path)
    donor_id = a.obs['Donor ID'].unique()[0]
    ADNC = a.obs["Overall AD neuropathological Change"].unique()[0]
    return a.obs["Subclass"].value_counts(), donor_id,ADNC

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_folder", type=str)
    parser.add_argument("output_folder", type=str)
    args = parser.parse_args()

    folder = args.input_folder
    file_paths = [os.path.join(folder, f) for f in os.listdir(folder) if f.endswith("RNA_multiome_subset.h5ad")]

    patient_status = dict()
    patient_subclass = dict()

    for patient_path in file_paths : 
        subclass,patient_ID,ADNC= analyze_patient(patient_path)
        patient_subclass[patient_ID] = subclass
        patient_status[patient_ID] = ADNC


    patient_subclass = pd.DataFrame(patient_subclass)
    patient_status  = pd.Series(patient_status)



    ADNC_COLORS = {
        "Not AD":       "#2196F3",  # bleu
        "Low":          "#4CAF50",  # vert
        "Intermediate": "#FF9800",  # orange
        "High":         "#F44336",  # rouge
    }



    col_colors = patient_status.map(ADNC_COLORS)
    patient_subclass.fillna(0, inplace=True)
    patient_subclass = patient_subclass.astype(int)
    patient_subclass_norm = patient_subclass.div(patient_subclass.sum(axis=0), axis=1) * 100

    g = sns.clustermap(patient_subclass_norm, 
                   col_colors=col_colors, 
                   annot=patient_subclass, 
                   fmt='d',
                   cmap="coolwarm",
                   figsize=(20, 12)) 

    #Adding the number of cell to the label
    reordered_cols = patient_subclass_norm.columns[g.dendrogram_col.reordered_ind]
    """
    g = sns.clustermap(patient_subclass, 
                   col_colors=col_colors, 
                   annot=patient_subclass, 
                   fmt='d',
                   cmap="coolwarm",
                   figsize=(20, 12)) 
    """
    #Adding the number of cell to the label
    #reordered_cols = patient_subclass.columns[g.dendrogram_col.reordered_ind]
    new_labels = [f"{col}\n(n={patient_subclass[col].sum()})" for col in reordered_cols]
    g.ax_heatmap.set_xticklabels(new_labels, rotation=45, ha='right', fontsize=7)

    legend_handles = [
        Line2D([0], [0], marker='s', color='w',
            markerfacecolor=color, markersize=11, label=label)
        for label, color in ADNC_COLORS.items()
    ]

    g.ax_col_dendrogram.legend(
    handles=legend_handles,
    title="ADNC",
    loc="center",
    frameon=False,
    )
    output_folder = args.output_folder
    os.makedirs(output_folder, exist_ok=True)

    g.savefig(os.path.join(output_folder, "clustermap.png"), dpi=150, bbox_inches="tight")