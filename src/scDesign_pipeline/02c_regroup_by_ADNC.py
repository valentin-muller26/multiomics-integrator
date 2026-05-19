import argparse
import anndata as ad
import gc
import pandas as pd
from pathlib import Path

"""
Groups donor-level h5ad files by ADNC neuropathological stage.

Usage:
    python grouping_by_ADNC.py <input_path> <output_path>
"""
parser = argparse.ArgumentParser()
parser.add_argument("input_path", type=str)
parser.add_argument("output_path", type=str)
# Optional argument to balance the groups by limiting each group to the size of the smallest non-empty group default is False
parser.add_argument(
        "--balance_groups",
        action="store_true",
        help="Limit each group to the size of the smallest non-empty group",
    )
args = parser.parse_args()


# Define the name of the metadata column containing the ADNC information
NAME_METADATA_ADNC = "Overall AD neuropathological Change"

# Dictionary linking the ADNC into group, allows for instance fusion of the ADNC in one group
ADNC_GROUP_MAP = {
    "Not AD":       "Not AD",
    "Low":          "Low",
    "Intermediate": "Intermediate",
    "High":         "High",
}

def grouping_by_ADNC(file_paths, assay_name, output_dir,balance_groups=False):
    """Groups and concatenates individual donor h5ad data files into ADNC group files based on ADNC neuropathological stage.

    Reads each AnnData file, retrieves its ADNC group from the obs metadata,
    concatenates all files belonging to the same group, and saves the merged
    result as a single .h5ad file per group. A manifest CSV listing the output
    paths is also written to the output directory.

    Args:
        file_paths (list[str]): List of paths to the input .h5ad files,
            one per donor.
        assay_name (str): Name of the assay (e.g. "RNA" or "ATAC"). Used
            as a subdirectory name under output_dir.
        output_dir (pathlib.Path): Root output directory where the grouped .h5ad
            files and the manifest CSV will be saved. Created if it does
            not exist.
        balance_groups (bool): If True, limits the number of files in each group to the size of the smallest non-empty group to ensure balanced group sizes. 
        Default is False.

    Returns:
        None.
"""

    # Store the paths of the files by group
    groups = {
        "Not AD":       [],
        "Low":          [],
        "Intermediate": [],
        "High":         [],
    }

    # Iterate over the files and retrieve the ADNC group for each file
    for patient_path in file_paths:
        print(f"Loading : {Path(patient_path).name}")
        adata = ad.read_h5ad(patient_path, backed="r")

        # Check if the ADNC group is present in the obs of the anndata object
        adnc_vals = adata.obs[NAME_METADATA_ADNC].unique()
        if len(adnc_vals) == 0:
            print(f"  Error: no ADNC values")
            adata.file.close()
            continue

        # Map the ADNC value to its group
        adnc_val = str(adnc_vals[0]).strip()
        group = ADNC_GROUP_MAP.get(adnc_val)

        if group is None:
            print(f"  WARN: Unknown ADNC value: '{adnc_val}', skipping")
            adata.file.close()
            continue

        print(f"  -> group : {group}  |  {adata.n_obs} cells")
        groups[group].append(patient_path)

        adata.file.close()
        gc.collect()

    manifest = {}
    
    # To ensure balanced groups, we will limit the number of files in each group to the size of the smallest group
    non_empty_counts = [len(donor_paths) for donor_paths in groups.values() if donor_paths]
    min_donor_count = min(non_empty_counts) if non_empty_counts else 0
    

    # Iterate over the groups, concatenate the files and save the result
    for group_name, donor_paths in groups.items():
        if not donor_paths:
            print(f"SKIP : No files in group '{group_name}'")
            continue
        
        # If balance_groups is True, limit the number of files in each group to the size of the smallest group
        if balance_groups:
            print(f"\n[{group_name}] {len(donor_paths)} donor(s) available")
            if len(donor_paths) > min_donor_count:
                donor_paths = donor_paths[:min_donor_count]
                groups[group_name] = donor_paths
                print(f"  -> Limiting to {min_donor_count} files for balance")

        # Load the data for concatenation
        adatas = [ad.read_h5ad(donor_path) for donor_path in donor_paths]
        merged = ad.concat(adatas, join="outer", merge="first", index_unique="-")
        del adatas
        gc.collect()

        # Add group metadata
        merged.obs["ADNC_group"] = group_name
        merged.uns["ADNC_group"] = group_name

        # Save the merged file for the group
        (output_dir / assay_name).mkdir(parents=True, exist_ok=True)
        h5ad_path = output_dir / assay_name / f"{group_name}.h5ad"
        merged.write_h5ad(h5ad_path)
        print(f"OK Saved -> {h5ad_path}  ({merged.n_obs} cells)")

        manifest[group_name] = h5ad_path
        del merged
        gc.collect()

    # Save the manifest CSV mapping each group name to its output path
    manifest_path = output_dir / assay_name / "manifest.csv"
    pd.Series(manifest).to_csv(manifest_path, index=True, header=False)

    # Summary of the groups created
    print("\n--- Summary ---")
    if not manifest:
        print("No groups were created.")
    else:
        print(f"Groups saved in {output_dir / assay_name}:")
        for group_name, h5ad_path in manifest.items():
            donor_list = [Path(p).name for p in groups[group_name]]
            print(f"  - {group_name} ({len(donor_list)} donors) -> {h5ad_path}")
            for donor in donor_list:
                print(f"      * {donor}")
        print(f"Manifest saved -> {manifest_path}")
    print("--- Done ---\n")


# Creating the output directory
output_dir = Path(args.output_path)
output_dir.mkdir(parents=True, exist_ok=True)

# Listing the donor RNA files
rna_paths = sorted(Path(args.input_path).glob("*RNA_multiome_subset.h5ad"))

# Listing the donor ATAC files
atac_paths = sorted(Path(args.input_path).glob("*ATAC_multiome_subset.h5ad"))

grouping_by_ADNC(rna_paths, "RNA", output_dir, balance_groups=args.balance_groups)
grouping_by_ADNC(atac_paths, "ATAC", output_dir, balance_groups=args.balance_groups)