# NOTE: code based on read_SA_seq.py -- parts that have been repurposed are annotated less
from collections import Counter
import io
import os
import re
import time
import pandas as pd
from Bio import Entrez, SeqIO

Entrez.email = "oscargarden27@gmail.com"


def find_specific_genes(genome_ids, target_genes=None, target_locus_tags=None, output_csv="SA_specific_genes.csv"):
    if target_genes is None:
        target_genes = []
    if target_locus_tags is None:
        target_locus_tags = []

    with open(output_csv, "w") as csv_file:
        csv_file.write(
            "Genome_ID,Genome_Name,Gene_Found,Sequence_Length,First_3,Last_3,GC_pct,AT_pct\n"
        )

        for genome_id in genome_ids:
            print(f"\nScanning genome {genome_id}...")

            try:
                # fetch the GenBank annotation file
                with Entrez.efetch(
                    db="nucleotide",
                    id=genome_id,
                    rettype="gb",
                    retmode="text",
                ) as fetch_handle:
                    gb_data = fetch_handle.read()

                gb_io = io.StringIO(gb_data)

                # parse through GenBank file
                for record in SeqIO.parse(gb_io, "genbank"):
                    short_desc = record.description[:40].replace(",", "")
                    print(f"  Successfully loaded: {short_desc}...")

                    found_any = False

                    # loop through features, trying to find CDS
                    for feature in record.features:
                        if feature.type == "CDS":
                            # extract identifiers separately
                            gene_names = feature.qualifiers.get("gene", [])
                            product_names = feature.qualifiers.get("product", [])
                            locus_tags = feature.qualifiers.get("locus_tag", [])

                            # create separate lowercase strings for independent lookups
                            names_str = " ".join(gene_names + product_names).lower()
                            loci_str = " ".join(locus_tags).lower()

                            match_found = False
                            matched_target = ""

                            # try to match by gene name first
                            for gene_target in target_genes:
                                if gene_target.lower() in names_str.split():
                                    match_found = True
                                    matched_target = gene_target
                                    break  # found match, stop looking through gene names

                            # try matching by locus tag
                            if not match_found:
                                for locus_target in target_locus_tags:
                                    if locus_target.lower() in loci_str.split():
                                        match_found = True
                                        matched_target = locus_target
                                        break  # found match, stop looking through locus tags

                            # extract and save data on success
                            if match_found:
                                found_any = True

                                # biopython auto extracts the exact sequence
                                gene_seq = str(feature.extract(record.seq)).upper()
                                total_bases = len(gene_seq)

                                # make sure we actually got a sequence
                                if total_bases == 0:
                                    continue

                                first_3 = gene_seq[:3]
                                last_3 = gene_seq[-3:]

                                # calc GC/AT percentages
                                counts = Counter(gene_seq)
                                gc_pct = ((counts["G"] + counts["C"]) / total_bases) * 100
                                at_pct = ((counts["A"] + counts["T"]) / total_bases) * 100

                                print(f"    [FOUND] {matched_target} | Length: {total_bases} | GC: {gc_pct:.2f}% | AT: {at_pct:.2f}%")

                                csv_file.write(f"{genome_id},{short_desc},{matched_target},{total_bases},{first_3},{last_3},{gc_pct:.2f},{at_pct:.2f}\n")

                    if not found_any:
                        print(
                            f"    [!] None of the target genes or locus tags were found in this genome."
                        )
            except Exception as e:
                print(f"  [Error] Failed to process ID {genome_id}: {e}")
            time.sleep(0.5)
    print(f"\nAnalysis complete! Results exported to '{output_csv}'")


# we need to separate the gene names and locus tags
df = pd.read_csv("Pseudomonas_sig_genes.csv")

# filter rows that are just a hyphen
active_df = df[df["gene_name"] != "-"]

# separate rows that are locus tags (for PA it would be PA followed by a number)
locus_mask = active_df["gene_name"].astype(str).str.contains(r"^PA\d+", case=False, na=False)

locus_list = active_df[locus_mask]["gene_name"].tolist()
genes_list = active_df[~locus_mask]["gene_name"].tolist()


find_specific_genes(
    genome_ids=["AE004091.2"],
    target_genes=genes_list,
    target_locus_tags=locus_list,
    output_csv="PAO1_specific_genes.csv"
)