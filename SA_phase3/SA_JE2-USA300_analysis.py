# NOTE: code based on SA_specific_genes.py -- parts that have been repurposed are annotated less

# NOTE: while coding, i ran into a GenBank vs RefSeq annotation rift/mismatch that I was able to work around:
# The 2006 original CP000255.1 (when researchers first sequenced USA300) was
# submitted to the GenBank database and the genes were named sequentially (SAUSA300_0001 through SAUSA300_2648).
#
# Then, much later when ncbi ran the og genome through their new automated pipeline (pgap), they were doing it to 
# standardize it for the RefSeq db. So, PGAP got rid of the 4 digit tags and replaced them w the 5 digit ones with "RS"
# (as a result, mapping the old SAUSA300_1759 to the new SAUSA300_RS09630). 
#
# This then requires the code to have NC_007793.1 (the RefSeq version) hardcoded at the bottom to read the RS tags in the .csv file


import csv
import io
import os
import re
import time
from Bio import Entrez, SeqIO
from collections import Counter

Entrez.email = "oscargarden27@gmail.com"

# reads first column of significant_list.csv automatically
# max_rows is the safety limiter during testing
def get_targets_from_csv(csv_filepath, max_rows=5):
    if not os.path.exists(csv_filepath):
        print(f"Error: could not locate '{csv_filepath}' in this directory.")
        return []
    
    targets = []
    # utf-8-sig strips the invisible byte order mark that excel puts on headers
    with open(csv_filepath, mode="r", encoding="utf-8-sig") as f:
        reader = csv.reader(f)
        header = next(reader, None)  # skip the header row (like "gene_name")

        for count, row in enumerate(reader):
            if max_rows and count >= max_rows:
                break
            # make sure the row isn't blank and grab col 1 (index 0)
            if row and row[0].strip():
                targets.append(row[0].strip())

    print(f"Ingested {len(targets)} target tags/symbols from CSV (Limited to {max_rows}).\n")
    return targets


def analyze_usa300_omni(genome_ids, input_csv, output_csv="SA_JE2-USA300_analysis.csv", row_limit=5):
    # first auto grab the targets (gene names in first col of sig list csv file)
    target_list = get_targets_from_csv(input_csv, max_rows=row_limit)
    if not target_list:
        return

    # Lowercase them matching purposes (case insensitive)
    targets_clean = [t.lower() for t in target_list]

    with open(output_csv, "w") as out_file:
        # added a 'Matched_Via' column so the CSV includes how it found it (like through gene symbol or locus tag etc)
        out_file.write("Genome_ID,Genome_Name,Target_Searched,Matched_Via,Sequence_Length,First_3,Last_3,GC_pct,AT_pct\n")

        for genome_id in genome_ids:
            print(f"Fetching GenBank map for {genome_id}...")

            try:
                # NOTE: IMPORTANT ADDITION: to force ncbi to give the actual dna seqs and gene feat.s 
                # for a RefSeq genome, rettype="gb" is changed to rettype="gbwithparts"
                with Entrez.efetch(db="nucleotide", id=genome_id, rettype="gbwithparts", retmode="text") as handle:
                    gb_data = handle.read()

                for record in SeqIO.parse(io.StringIO(gb_data), "genbank"):
                    short_desc = record.description[:40].replace(",", "")
                    print(f"  Loaded: {short_desc}...")

                    # keep track of findings to report what was missing at the end
                    found_targets = set()

                    for feature in record.features:
                        # look at 'CDS' too !! (its where the real dna sequence lives...)
                        if feature.type in ["CDS", "gene", "rRNA", "tRNA"]:

                            # open all three possibilities ncbi keeps identifiers in
                            g_names   = [x.lower() for x in feature.qualifiers.get("gene", [])]
                            l_tags    = [x.lower() for x in feature.qualifiers.get("locus_tag", [])]
                            old_tags  = [x.lower() for x in feature.qualifiers.get("old_locus_tag", [])]

                            # Combine the locus tags into one 'pool'
                            all_locus_pool = l_tags + old_tags

                            for target_orig, target_low in zip(target_list, targets_clean):
                                if target_orig in found_targets:
                                    continue  # this means the seq was already found for this item

                                match_type = None

                                # CHECK 1: Did it match a gene symbol? (like: 'mtlA')
                                if target_low in g_names:
                                    match_type = f"Gene Symbol ({feature.qualifiers['gene'][0]})"

                                # CHECK 2: Did it match an exact locus tag? ('SAUSA300_RS09630')
                                elif target_low in all_locus_pool:
                                    match_type = "Exact Locus Tag"

                                # CHECK 3: 'catch-all' substring check 
                                # (this basically catches if the .csv file has 'RS09630', but genbank holds 'SAUSA300_RS09630')
                                else:
                                    for tag in all_locus_pool:
                                        if target_low in tag:
                                            match_type = f"Sub-string Locus ({tag.upper()})"
                                            break

                                # CHECK 4: pure number matching (aka how to bypass "RS" issues in the locus tag)
                                # this catches if the csv asks for (for example) RS09630 but genbank says '09630'
                                if not match_type and "_" in target_low:
                                    # strip all letters from teh right side of the underscore, leaving only digits!
                                    target_numbers = re.sub(r'\D', '', target_low.split("_")[-1])
                                    
                                    if target_numbers: # Make sure there are actually numbers to compare
                                        #loops through genbanks locus tags
                                        for tag in all_locus_pool: 
                                            if "_" in tag:
                                                tag_numbers = re.sub(r'\D', '', tag.split("_")[-1])

                                                if target_numbers == tag_numbers:
                                                    match_type = f"Pure Number Match ({tag.upper()})"
                                                    break 

                                # if any of the 4 checks passed:
                                if match_type:
                                    found_targets.add(target_orig)

                                    gene_seq = str(feature.extract(record.seq)).upper()
                                    length = len(gene_seq)

                                    if length == 0:
                                        continue

                                    first_3 = gene_seq[:3]
                                    last_3 = gene_seq[-3:]
                                    counts = Counter(gene_seq)
                                    gc = ((counts['G'] + counts['C']) / length) * 100
                                    at = ((counts['A'] + counts['T']) / length) * 100

                                    print(f"    [FOUND] {target_orig} via {match_type} | Len: {length}bp | GC: {gc:.1f}% | AT: {at:.1f}% | {first_3}...{last_3}")

                                    out_file.write(f"{genome_id},{short_desc},{target_orig},{match_type},{length},{first_3},{last_3},{gc:.2f},{at:.2f}\n")
                                    break

                    # report what slipped through the cracks !!
                    missing = set(target_list) - found_targets
                    if missing:
                        print(f"    [!] Failed to find sequences for: {', '.join(missing)}")

            except Exception as e:
                print(f"  [Error processing {genome_id}]: {e}")

            time.sleep(0.5)

    print(f"\nDone! Output saved to '{output_csv}'")

# Run the script:
# can pass a single ID or a whole list of IDs here -- 
# since this is looking specifically at Staphylococcus aureus strain JE2 USA300, its GenBank ID is CP092052.1
# so that is hardcoded in for test_genome_ids

usa300_accession = ["NC_007793.1"]

my_csv_file = "significant_list.csv" # input csv file (for first row reading automation)

# Set row_limit=5 to test. When ready for the whole spreadsheet, change to: row_limit=None
analyze_usa300_omni(usa300_accession, my_csv_file, output_csv="USA300_Results.csv", row_limit=None)

# find_specific_genes(test_genome_ids, target_genes_list)