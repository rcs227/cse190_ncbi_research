import io
import time
from Bio import Entrez, SeqIO
from collections import Counter

# Note: you need to tell NCBI who you are because they may block you if youre making lots
# of requests at once. Bio.Entrez associates this email with each call to Entrez.
Entrez.email = "oscargarden27@gmail.com"
# Note: suggested to get API key, as this project handles large amounts of data.
# You can make at most 10 queries per second with it; otherwise, the cap is at 3.
# Entrez.api_ley = "MyAPIkey"

# attempts to search ncbi, retrying if the backend fails.
def ncbi_search(search_term, max_results, retries=3, delay=2):
    # need to include retries and delay because of Search Backend Failed 
    # runtime error (has to do with too many queries at once to ncbi).

    for attempt in range(retries):
        try:
            with Entrez.esearch(db="nucleotide", term=search_term, retmax=max_results) as handle:
                return Entrez.read(handle)
        except RuntimeError as e:
            if "Search Backend failed" in str(e) and attempt < retries - 1:
                print(f"   [NCBI Server Hiccup] Backend failed. Attempt {attempt + 1}/{retries}. Retrying in {delay} seconds...")
                time.sleep(delay)
                delay *= 2  #exponential backoff
            else:
                raise e

# searches for S. aureus genomes, streams them, and calculates nucleotide %
# without saving the big sequence files to disk.
def stream_and_analyze_staph_genomes(max_results=5):
    #Note: max_results=5 is testing number to see what code and db can handle

    print("Searching NCBI for Staphylococcus aureus genomes...")

    # first, search the Nucleotide database for S. aureus complete genomes
    search_term = "txid1280[Organism] AND complete genome[Title]"
    #Note: use the exact NCBI Taxonomy ID for S. aureus (txid1280)
    # db backends prefer numbers since it makes the query faster and less likely to fail

    try:
        search_results = ncbi_search(search_term, max_results)
        id_list = search_results["IdList"]
    except RuntimeError as e:
        print(f"\nFailed to connect to NCBI: {e}")
        return

    print(f"Found {len(id_list)} matches. Beginning streamed analysis...\n")

    # Open a summary file to save results
    with open("staph_genome_summary.csv", "w") as summary_file:
        summary_file.write("Accession,Description,Total_Bases,A_pct,T_pct,G_pct,C_pct,GC_pct,AT_pct\n")

        #next, loop through all genome IDs to stream the data
        for i, genome_id in enumerate(id_list):
            print(f"[{i+1}/{len(id_list)}] Streaming and analyzing genome ID: {genome_id}...")

            try:
                # note: efetch requests the specific genome in fasta format
                with Entrez.efetch(db="nucleotide", id=genome_id, rettype="fasta", retmode="text") as fetch_handle:
                    #read the web response into memory as text
                    fasta_data = fetch_handle.read()
                
                # Convert the raw text into an in-memory obj so SeqIO can read it
                fasta_io = io.StringIO(fasta_data)

                for record in SeqIO.parse(fasta_io, "fasta"):
                    # need to convert sequence to uppercase string in case of inconsistencies
                    sequence_str = str(record.seq).upper()
                    # counter counts all characters at once!
                    counts = Counter(sequence_str)

                    # total strictly of A, T, G, C
                    total_atgc = counts['A'] + counts['T'] + counts['G'] + counts['C']
                    if total_atgc == 0:
                        continue
                    
                    # Calculate percentages
                    a_pct = (counts['A'] / total_atgc) * 100
                    t_pct = (counts['T'] / total_atgc) * 100
                    g_pct = (counts['G'] / total_atgc) * 100
                    c_pct = (counts['C'] / total_atgc) * 100
                    gc_content = g_pct + c_pct
                    at_content = a_pct + t_pct
                    
                    short_desc = record.description[:40].replace(",", "") 

                    # write to local CSV file !!
                    summary_file.write(f"{record.id},{short_desc},{total_atgc},{a_pct:.2f},{t_pct:.2f},{g_pct:.2f},{c_pct:.2f},{gc_content:.2f},{at_content:.2f}\n")
                    
                    #note, only prints GC and AT content for each, go to csv file for separate percentages
                    print(f"    Success: {short_desc}... | GC Content: {gc_content:.2f}%")
                    print(f"    Success: {short_desc}... | AT Content: {at_content:.2f}%")
            except Exception as e:
                print(f"    Error processing ID {genome_id}: {e}")
            
            # pause between requests to prevent server strain -- NCBI rules
            time.sleep(0.5) 

    print("\nAnalysis complete! Results saved to 'staph_genome_summary.csv'")

# Run the function (configured to fetch the first 5 for testing)
stream_and_analyze_staph_genomes(max_results=5)

                

    
    