import pandas as pd
import os
import bvbrc as bv

client = bv.GenomeClient()
# extract s. aureus accessions
sa_df = pd.read_csv("SA_local_ids_summary.csv")
sa_accessions = sa_df["Accession"].tolist()
df = pd.DataFrame()
filename = "staph_from_accessions.csv"
for accession in sa_accessions:
    file_exists = os.path.exists(filename)

    # check refseq accession
    q = bv.query(
        species = "Staphylococcus aureus",
        refseq_accessions = accession
    )
    response = client.submit_query(q)
    df = response.to_pandas()

    # if empty, check genbank accession
    if df.empty:
        q = bv.query(
            species = "Staphylococcus aureus",
            genbank_accessions = accession
        )
        response = client.submit_query(q)
        df = response.to_pandas()
    else:
        print(f":) Found data for genome with RefSeq accession {accession}\n")

    # check if GenBank was found, if not then the data is not available
    if df.empty:
        print(f"    :( No data found for genome with RefSeq/GenBank accession {accession}\n")
        continue
    else:
        print(f":) Found data for genome with GenBank accession {accession}\n")

    df.to_csv(filename, mode="a", header=not file_exists, index=False)