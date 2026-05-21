"""
Analysis script for MapMyCells mapping results.

Loads a MapMyCells output CSV, computes classification accuracy (global and
per cell type), and generates confusion matrices (normalized and raw).

Usage:
    python analysis_mapmycells.py <file_name> <input_path> <output_path>
"""
#mapmycells
import pandas as pd
import argparse
from pathlib import Path
from sklearn import metrics
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser(description="Analyze MapMyCells output and generate classification statistics and confusion matrices.")
parser.add_argument("file_name", type=str, help="Name of the input files.")
parser.add_argument("input_path", type=str, help="Path to the input CSV file containing MapMyCells results.")
parser.add_argument("output_path", type=str, help="Path to the output directory for the analysis results.")
args = parser.parse_args()

def plot_confusion_matrix(prediction, labels, output_path, normalize=None):
    """Generate and save a confusion matrix plot from prediction data.

    Args:
        prediction: DataFrame with columns 'cell_type' (truth) and 'subclass_name' (prediction).
        labels: List of class labels (sorted) for consistent axis ordering.
        output_path: Path where the .png will be saved.
        normalize: Passed to sklearn's confusion_matrix. Use "true" for row-normalized
                   (percentages) or None for raw counts.
    """
    cm = metrics.confusion_matrix(
        prediction["cell_type"], prediction["subclass_name"],
        labels=labels, normalize=normalize
    )
    cm_display = metrics.ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=labels)
    fig, ax = plt.subplots(figsize=(14, 12))
    cm_display.plot(ax=ax, cmap="Blues", xticks_rotation=90)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


#Create the output directory 
output_dir = Path(args.output_path)
output_dir.mkdir(parents=True, exist_ok=True)

#Load the results of mapmycells and convert the cell_id to the cell type by taking the first part of the cell_id before the "_"
# The true cell type is the substring before the first underscore
# (e.g. "Astrocyte_001" -> "Astrocyte").
mapmycells_data = pd.read_csv(args.input_path, skiprows=4)  # skip the first 4 rows of the csv file which contain metadata
mapmycells_data["cell_type"] = mapmycells_data["cell_id"].str.split("_").str[0]
number_cell = mapmycells_data.shape[0]

# Global classification accuracy
percentage_classification = (sum(mapmycells_data["cell_type"] == mapmycells_data["subclass_name"]) / number_cell) * 100
print("Percentage correct classification: ", percentage_classification)

# Save the global accuracy to a text file
global_summary_path = output_dir / f"{args.file_name}_global_accuracy.csv"
pd.DataFrame([{
    "file_name": args.file_name,
    "n_cells": number_cell,
    "global_accuracy_percent": percentage_classification
}]).to_csv(global_summary_path, index=False)
print(f"Global accuracy saved to: {global_summary_path}")

# Per cell type classification %
percentage_per_celltype = mapmycells_data.groupby("cell_type").apply(
    lambda g: (g["cell_type"] == g["subclass_name"]).mean() * 100
)
percentage_per_celltype.to_csv(output_dir / f"{args.file_name}_percentage_per_celltype.csv", header=["percentage_correct"])

# Confusion matrices (percentage and raw)
prediction = mapmycells_data[["cell_type", "subclass_name"]]
labels = sorted(mapmycells_data["cell_type"].unique())

plot_confusion_matrix(
    prediction, labels,
    output_dir / f"{args.file_name}_confusionmatrix_percentage.png",
    normalize="true"
)
plot_confusion_matrix(
    prediction, labels,
    output_dir / f"{args.file_name}_confusionmatrix_raw.png",
    normalize=None
)