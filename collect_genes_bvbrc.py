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
        print(f"Found {len(gene_df)} antibiotic resistance genes for [{id}] {name}.")

        # filter for antibiotic resistance property
        if "property" in gene_df.columns:
            gene_df = gene_df[
                gene_df["property"].str.contains(
                    "Antibiotic Resistance", case=False, na=False
                )
            ]

        if gene_df.empty:
            print(f"    No antibiotic resistance genes found for [{id}] {name}.")
            continue

        gene_df.to_csv(gene_file_path, mode="a", index=False, header=not os.path.exists(gene_file_path))
        print(f":) Found {len(gene_df)} antibiotic resistance genes for [{id}] {name}!")

    print(f"All genomes searched and saved to {gene_file_path}!")


collect_genes("Staphylococcus aureus")