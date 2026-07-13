import pandas as pd

import csv
import re

def extract_shine_dalgarno_regions(dna_filepath, genes_list, output_csv, sd_sequence="AGGAGG"):
    """
    Reads a single-line DNA sequence, finds specified genes, and looks for 
    Shine-Dalgarno sequences in the 100 bases upstream of each gene.
    Saves the extracted regions (15 bases + SD + 15 bases) to a CSV.
    """
    
    # 1. Read the DNA sequence from the text file
    with open(dna_filepath, 'r') as file:
        # .strip() and .replace() ensure we remove any hidden newlines/spaces
        full_dna = file.read().strip().replace('\n', '')
        
    dna_length = len(full_dna)

    # 2. Open the CSV file for writing
    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        # Write the header
        writer.writerow(["Gene_Name", "Extracted_SD_Region"])

        # 3. Iterate through each gene in the provided list
        for gene in genes_list:
            gene_name = gene['name']
            gene_seq = gene['sequence']
            
            # Find all starting positions of the gene in the full DNA sequence
            # (Allows for multiple copies of the same gene)
            for gene_match in re.finditer(gene_seq, full_dna):
                gene_start = gene_match.start()
                
                # Define the 100-character window before the gene
                # max(0, ...) ensures we don't get a negative index if the gene is near the start
                upstream_start = max(0, gene_start - 100)
                upstream_region = full_dna[upstream_start:gene_start]
                
                # 4. Search for the Shine-Dalgarno sequence within this 100-base window
                for sd_match in re.finditer(sd_sequence, upstream_region):
                    # Calculate the global coordinates of the SD sequence in the full DNA string
                    sd_global_start = upstream_start + sd_match.start()
                    sd_global_end = upstream_start + sd_match.end()
                    
                    # 5. Extract 15 bases before and 15 bases after (including the SD itself)
                    # max() and min() prevent index out-of-bounds errors at the ends of the genome
                    extract_start = max(0, sd_global_start - 15)
                    extract_end = min(dna_length, sd_global_end + 15)
                    
                    extracted_string = full_dna[extract_start:extract_end]
                    
                    # Save to CSV
                    writer.writerow([gene_name, extracted_string])


def find_number_genes():
    csv_path = 'PAO1_DNA.csv'

    df = pd.read_csv(csv_path)
    rows = df.itertuples(index=False)

    # 2. Read the text file and remove newlines to make it one continuous string
    text_path = 'pao1_sequence.txt'
    with open(text_path, 'r', encoding='utf-8') as file:
        # .read() gets the whole file, .replace('\n', '') glues the lines together
        # We also replace '\r' just in case it's a Windows-formatted text file
        target_text = file.read()

    # 3. Search for each string and count the matches
    total_expected = len(df)

    genes = []

    for row in rows:
        item = row.Full_Sequence
        name = row.gene_name
        if item in target_text:
            genes.append({"name": name, "sequence": item})

    # 4. Print the final results
    print(f"Match Results: Found {len(genes)} out of {total_expected} expected strings.")
    return genes

def remove_newlines_from_file(file_path, output_path=None):
    """
    Reads a text file, removes all newline characters, and writes it back out.
    """
    # 1. Read the entire content into memory first
    with open(file_path, 'r', encoding='utf-8', errors='ignore') as file:
        content = file.read()
    
    # 2. Strip out the newline characters
    clean_content = content.replace('\r', '').replace('\n', '')
    
    # 3. Figure out the destination
    target_path = output_path if output_path else file_path
    
    # 4. Write out the clean content only after the read is fully finished
    with open(target_path, 'w', encoding='utf-8') as file:
        file.write(clean_content)
        
    print(f"Successfully formatted! Saved to: {target_path}")

genes = find_number_genes()


extract_shine_dalgarno_regions(dna_filepath="pao1_sequence.txt", genes_list=genes, output_csv="shine_res.csv", sd_sequence="GGAG|GAGG")