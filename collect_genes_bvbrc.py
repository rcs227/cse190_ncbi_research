import pandas as pd
import bvbrc as bv
import os

def collect_genes(species: str):
    # path to file that contains every name of genome of the species
    id_file_path = "bvbrc_" + species.replace(" ", "_") + "_genomes.csv"
    id_file_exists = os.path.exists(id_file_path)

    # if we don't have the IDs yet collect them
    if not id_file_exists:
        print("Collecting genome names and IDs from BV-BRC database...")
        client = bv.GenomeClient()
        q = bv.query(species = species, limit="max")
        response = client.submit_query(q)
        df = response.to_pandas()
        if df.empty:
            print("     No genomes found.")
        else:
            print(f"Found {len(df)} genomes for {species}.")
        df[["genome_id", "genome_name"]].to_csv(id_file_path, index=False)
        print(f"Saved genomes to {id_file_path}.")
    
    # get list of genomes
    id_df = pd.read_csv(id_file_path, dtype={"genome_id": str})
    client = bv.SpGeneClient()
    gene_file_path = species.replace(" ", "_") + "_gene_info.csv"

    for id, name in zip(id_df['genome_id'], id_df['genome_name']):
        q = bv.query(genome_id = id)

        # retrieve all specialty genes of genome
        print(f"Querying specialty genes for [{id}] {name}...")
        response = client.submit_query(q)
        gene_df = response.to_pandas()

        if gene_df.empty:
            print(f"     No specialty genes found for [{id}] {name}.\n")
            continue
        print(f"Found {len(gene_df)} specialty genes for [{id}] {name}.")

        # filter for antibiotic resistance property
        if "property" not in gene_df.columns:
            print(f"    No 'property' column found for [{id}] {name}, skipping.")
            continue

        gene_df = gene_df[gene_df["property"].str.contains("antibiotic resistance", case=False, na=True)]

        if gene_df.empty:
            print(f"    No antibiotic resistance genes found for [{id}] {name}.")
            continue

        if os.path.exists(gene_file_path):
            # we can't just append each row as we go without checking the data because some column data is left blank.
            # this causes a shift in the data. for example, if "date added" is left blank, 
            # the next column would fill in its place even though it would be incorrect.

            # read existing headers and merge with new columns
            existing_cols = pd.read_csv(gene_file_path, nrows=0).columns.tolist()
            all_cols = existing_cols + [c for c in gene_df.columns if c not in existing_cols]

            # rewrite the file with any new columns added, filling old rows with NaN
            if all_cols != existing_cols:
                existing_df = pd.read_csv(gene_file_path)
                existing_df = existing_df.reindex(columns=all_cols)
                existing_df.to_csv(gene_file_path, index=False)

            # align new rows to the full column set before appending
            gene_df = gene_df.reindex(columns=all_cols)
            gene_df.to_csv(gene_file_path, mode="a", index=False, header=False)
        else:
            # first write — just dump as-is
            gene_df.to_csv(gene_file_path, index=False)
        print(f":) Found {len(gene_df)} antibiotic resistance genes for [{id}] {name}!\n")

    print(f"All genomes searched and saved to {gene_file_path}!")


collect_genes("Staphylococcus aureus")