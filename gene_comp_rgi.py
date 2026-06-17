import json

def calculate_dna_stats(sequence):
    # Normalize to uppercase to ensure accurate counting
    seq = sequence.upper()
    total_len = len(seq)
    if total_len == 0:
        return 0, 0, "", ""
    
    # Calculate counts
    g_count = seq.count('G')
    c_count = seq.count('C')
    a_count = seq.count('A')
    t_count = seq.count('T')
    
    gc_percent = ((g_count + c_count) / total_len) * 100
    at_percent = ((a_count + t_count) / total_len) * 100
    
    # Extract first and last 3 bases
    first_3 = seq[:3]
    last_3 = seq[-3:]
    
    return round(at_percent, 2), round(gc_percent, 2), first_3, last_3

def process_genome_data(input_json):
    organized_data = {}

    for key, value in input_json.items():
        # Extract accession from the key (e.g., "CP180884.1_...")
        # Splitting by '_' and taking the first part
        accession = key.split('_')[0]
        
        if accession not in organized_data:
            organized_data[accession] = {}
            
        for gene_id, gene_data in value.items():
            dna = gene_data.get("dna_sequence_from_broadstreet", "")
            at_pct, gc_pct, first, last = calculate_dna_stats(dna)
            
            organized_data[accession][gene_id] = {
                "name": gene_data.get("model_name"),
                "at_percent": at_pct,
                "gc_percent": gc_pct,
                "sequence_ends": f"{first} {last}"
            }
            
    return organized_data

# Usage
# with open('input.json', 'r') as f:
#     data = json.load(f)
# result = process_genome_data(data)
# print(json.dumps(result, indent=4))

with open("test.json", "r") as f:
    data = json.load(f)
result = process_genome_data(data)
print(json.dumps(result, indent=4))