import argparse
import anndata as ad
import gc
import pandas as pd
from pathlib import Path
import os

#Dictionnary linking the ADNC into group allow for instance fusion of the ADNC in on group
ADNC_GROUP_MAP = {
    "Not AD":       "Not AD",
    "Low":          "Low",
    "Intermediate": "Intermediate",
    "High":         "High",
}

def grouping_by_ADNC(file_paths,assay_name,output_dir):
    # store the path of the file by groupe
    groups = {
        "Not AD":       [],
        "Low":          [],
        "Intermediate": [],
        "High":         [],
    }

    for patient_path in file_paths:
        print(f"Loading : {Path(patient_path).name}")
        # Opening in backed mod for reading the ADNC group
        a = ad.read_h5ad(patient_path, backed="r")

        adnc_vals = a.obs["Overall AD neuropathological Change"].unique()
        if len(adnc_vals) == 0:
            print(f"  [Error] no ADNC values")
            a.file.close()
            continue

        adnc_val = str(adnc_vals[0]).strip()
        group = ADNC_GROUP_MAP.get(adnc_val)

        if group is None:
            print(f"  [WARN] Unknown ADNC value: '{adnc_val}', skipping")
            a.file.close()
            continue

        print(f"  -> group : {group}  |  {a.n_obs} cells")
        groups[group].append(patient_path)  # add the path to the file
        a.file.close()
        gc.collect()

    list_files_saved = {}
    # Fuse group by group
    for group_name, paths in groups.items():
        if not paths:
            print(f"[SKIP] No files in group '{group_name}'")
            continue

        print(f"Fusing {len(paths)} file(s) for '{group_name}'...")

        # loading the data only for the concatenation 
        adatas = [ad.read_h5ad(p) for p in paths]
        merged = ad.concat(adatas, join="outer", merge="first", index_unique="-")
        del adatas
        gc.collect()

        merged.obs["ADNC_group"] = group_name
        merged.uns["ADNC_group"] = group_name

        (output_dir / assay_name).mkdir(parents=True, exist_ok=True)
        out_path = output_dir / assay_name / f"{group_name}.h5ad"
        merged.write_h5ad(out_path)
        print(f"[OK] Saved -> {out_path}  ({merged.n_obs} cells)")
        #list of group created with the path
        list_files_saved[group_name] = out_path
        del merged
        gc.collect()
    
    out_path = output_dir / assay_name / "manifest.csv"
    list_files_saved = pd.Series(list_files_saved)
    list_files_saved.to_csv(out_path, index=True,header=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_path", type=str)
    parser.add_argument("output_path", type=str)
    args = parser.parse_args()

    #Creating the output directory
    output_dir = Path(args.output_path)
    output_dir.mkdir(parents=True, exist_ok=True)

    #Listing the donors RNA files
    rna_path = [
        os.path.join(args.input_path, f)
        for f in os.listdir(args.input_path)
        if f.endswith("RNA_multiome_subset.h5ad")
    ]
    #Listing the donors ATAC files
    atac_path = [
        os.path.join(args.input_path, f)
        for f in os.listdir(args.input_path)
        if f.endswith("ATAC_multiome_subset.h5ad")
    ]
    grouping_by_ADNC(rna_path,"RNA",output_dir)
    grouping_by_ADNC(atac_path,"ATAC",output_dir)


    print("Finished")