import anndata
import mygene
import argparse
from pathlib import Path
from os.path import basename
import os 

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_rna_path", type=str)
    parser.add_argument("output_path", type=str)
    args = parser.parse_args()

    #Create the output directory 
    output_dir = Path(args.output_path)
    output_dir.mkdir(parents=True, exist_ok=True)

    #Extract the file name and extension
    name = basename(args.input_rna_path)
    fileName, fileExtension = os.path.splitext(name)

    #Load the h5ad
    a = anndata.read_h5ad(args.input_rna_path)
    mg = mygene.MyGeneInfo()

    #query to convert the gene name to ensembl ID
    results = mg.querymany(
        a.var_names.tolist(),
        scopes='symbol,alias,other_names',
        fields='ensembl.gene,symbol',
        species='human',
        returnall=True,
        as_dataframe=True
    )


    #Extract the mapping from the results and apply it to the var_names of the anndata object
    df = results['out'].copy()

    df['ensembl_id'] = df['ensembl.gene'].apply(
        lambda x: x if isinstance(x, str) else (x[0] if isinstance(x, list) else None)
    )
    mapping = df[df['ensembl_id'].notna()].groupby('query')['ensembl_id'].first()
    a.var['ensembl_id'] = a.var_names.map(mapping)

    #Show statistic 
    print(f"Mapped: {a.var['ensembl_id'].notna().sum()} / {len(a.var_names)}")
    print(f"Unmapped: {a.var['ensembl_id'].isna().sum()}")

    #Replace the gene name by the ensembl ID and save the new h5ad file
    a_mapped = a[:, a.var['ensembl_id'].notna()].copy()
    a_mapped.var_names = a_mapped.var['ensembl_id'].values
    a_mapped.var_names_make_unique()

    output_file = output_dir / f"{fileName}_ens{fileExtension}"
    #write the newfile
    a_mapped.write_h5ad(output_file)
