# NOTE: code based on read_SA_seq.py -- parts that have been repurposed are annotated less
import io
import os
import re
import time
from Bio import Entrez, SeqIO
from collections import Counter

Entrez.email = "oscargarden27@gmail.com"

def find_specific_genes(genome_ids, target_genes, output_csv="SA_specific_genes.csv"):
    #convert target list to lowercase for matching purposes
    target_genes_lower = [g.lower() for g in target_genes]

    #first open the csv file to write the exported results
    with open(output_csv, "w") as csv_file:
        # write the header row
        csv_file.write("Genome_ID,Genome_Name,Gene_Found,Sequence_Length,First_3,Last_3,GC_pct,AT_pct\n")

        for genome_id in genome_ids:
            print(f"\nScanning genome {genome_id} for {target_genes}...")
        
        try:
            # NOTE: the rettype should be "gb" (aka GenBank) instead of "fasta"
            # doesnt download the raw letters, instead downloads annotated map of genome
            # here, were fetching the GenBank annotation file
            with Entrez.efetch(db="nucleotide", id=genome_id, rettype="gb", retmode="text") as fetch_handle:
                gb_data = fetch_handle.read()
                
            gb_io = io.StringIO(gb_data)
            
            # parse through GenBank file
            for record in SeqIO.parse(gb_io, "genbank"):
                short_desc = record.description[:40].replace(",", "")
                print(f"  Successfully loaded: {short_desc}...")
                
                found_any = False
                
                #check BOTH CDS and gene features to avoid ncbi annotation quirks
                # Loop through all the "features" (aka annotated parts of the genome)
                for feature in record.features:
                    if feature.type in ["CDS", "gene"]:
                        # extract gene name or product name
                        gene_names = feature.qualifiers.get("gene", [])
                        product_names = feature.qualifiers.get("product", [])
                        
                        #combine them into one string to search through
                        all_names = " ".join(gene_names + product_names).lower()

                        # check if any of the target genes are explicitly in the name/product
                        for target in target_genes_lower:
                            # NOTE: use exact word matching to ensure stuff like "rpoC-like" isnt accidentally grabbed
                            if target in all_names.split():
                                found_any = True
                            
                                # biopython auto extracts the exact sequence for the specific gene
                                gene_seq = str(feature.extract(record.seq)).upper()
                                total_bases = len(gene_seq)

                                # skip if something went wrong and sequence is empty
                                if total_bases == 0:
                                    continue
                                
                                #get first and last 3 bases
                                first_3 = gene_seq[:3]
                                last_3 = gene_seq[-3:]
                            
                                # calc GC/AT exactly like it was done before
                                counts = Counter(gene_seq)
                                gc_pct = ((counts['G'] + counts['C']) / total_bases) * 100
                                at_pct = ((counts['A'] + counts['T']) / total_bases) * 100
                                                            
                                # print to terminal
                                print(f"    [FOUND] {target.upper()} | Length: {total_bases} | GC: {gc_pct:.2f}% | AT: {at_pct:.2f}% | Bases: {first_3}...{last_3}")
                                    
                                # write to csvv
                                csv_file.write(f"{genome_id},{short_desc},{target.upper()},{total_bases},{first_3},{last_3},{gc_pct:.2f},{at_pct:.2f}\n")
                                    
                                # once the feature is found, break the target loop to avoid double-counting
                                break
                
                if not found_any:
                    print(f"    [!] None of the target genes were found under standard names in this genome.")
        except Exception as e:
            print(f"  [Error] Failed to process ID {genome_id}: {e}")
        time.sleep(0.5) # pause for ncbi servers
    print(f"\nAnalysis complete! Results exported to '{output_csv}'")


# --- Run the Script ---
# can pass a single ID or a whole list of IDs here
# TODO: extract the ids from SA_local_ids_summary.csv 
test_genome_ids = ["CP180884.1"] 
target_genes_list = ["rpoC", "adk", "tdk", "apt", "hpt", "pstC"]

find_specific_genes(test_genome_ids, target_genes_list)