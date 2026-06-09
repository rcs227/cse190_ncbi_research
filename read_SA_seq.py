import io
import os
import time
from Bio import Entrez, SeqIO
from collections import Counter

# Note: you need to tell NCBI who you are because they may block you if youre making lots
# of requests at once. Bio.Entrez associates this email with each call to Entrez.
Entrez.email = "oscargarden27@gmail.com"
# Note: suggested to get API key, as this project handles large amounts of data.
# You can make at most 10 queries per second with it; otherwise, the cap is at 3.
# Entrez.api_ley = "MyAPIkey"

# reads NCBI genome ids from local txt file and streams their data. 
#'max-to-process' allows there to be a limit to the run (like testing w 5 ids first)
def stream_SA_local_id_list(id_file_path, max_to_process=None):
    if not os.path.exists(id_file_path):
        print(f"Error: The file '{id_file_path}' was not found in this folder")
        return
    
    # first, read ids from local txt file
    print(f"Reading IDs from {id_file_path}...")
    with open(id_file_path, "r") as f:
        #read lines, strip whitespace/newlines and ignore empty lines if there are any
        id_list = [line.strip() for line in f if line.strip()]
    
    total_ids = len(id_list)
    print(f"Loaded {total_ids} IDs from file.")

    #can slice the list if want to test a few first
    if max_to_process:
        id_list = id_list[:max_to_process]
        print(f"Limiting execution to the first {len(id_list)} IDs for this run.\n")
    else:
        print("Beginning analysis on all IDs...\n")

    #next, open the summary csv to save results on the fly
    output_csv = "SA_local_ids_summary.csv"
    with open(output_csv, "w") as summary_file:
        summary_file.write("Accession,Description,Total_Bases,A_pct,T_pct,G_pct,C_pct,GC_pct,AT_pct\n")

        #next, loop through the local ids and stream from ncbi
        for i, genome_id in enumerate(id_list):
            #next, loop through all genome IDs to stream the data
            print(f"[{i+1}/{len(id_list)}] Streaming and analyzing genome ID: {genome_id}...")

            try:
                # fetch the specific genome structure using the current id
                with Entrez.efetch(db="nucleotide", id=genome_id, rettype="fasta", retmode="text") as fetch_handle:
                    fasta_data = fetch_handle.read()
                
                #parse the downloaded text entirely in mem
                fasta_io = io.StringIO(fasta_data)

                for record in SeqIO.parse(fasta_io, "fasta"):
                    sequence_str = str(record.seq).upper()
                    counts = Counter(sequence_str)

                    total_atgc = counts['A'] + counts['T'] + counts['G'] + counts['C']
                    if total_atgc == 0:
                        print(f"    Warning: No valid sequences found for ID {genome_id}")
                        continue
                        
                    # Calculate percentages
                    a_pct = (counts['A'] / total_atgc) * 100
                    t_pct = (counts['T'] / total_atgc) * 100
                    g_pct = (counts['G'] / total_atgc) * 100
                    c_pct = (counts['C'] / total_atgc) * 100
                    gc_content = g_pct + c_pct
                    at_content = a_pct + t_pct
                    
                    # clean up description text for csv safety
                    short_desc = record.description[:50].replace(",", "")
                    
                    # Write row to csv
                    summary_file.write(f"{record.id},{short_desc},{total_atgc},{a_pct:.2f},{t_pct:.2f},{g_pct:.2f},{c_pct:.2f},{gc_content:.2f},{at_content:.2f}\n")
                    print(f"    Success: {short_desc}... | GC Content: {gc_content:.2f}%")
                    print(f"    Success: {short_desc}... | AT Content: {at_content:.2f}%")
            
            except Exception as e:
                # If NCBI hiccups on one specific ID, log it and keep going
                print(f"    [Error] Failed to process ID {genome_id}: {e}")
            
            # Brief pause for NCBI guidelines
            time.sleep(0.5) 

    print(f"\nAnalysis complete! Check your folder for '{output_csv}'")

# Set max_to_process=None when ready to let it run through the whole file!
stream_SA_local_id_list("SA_nucleotide_genomes.txt", max_to_process=10)
