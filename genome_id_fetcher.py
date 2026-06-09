from Bio import Entrez

Entrez.email = "rcs227@lehigh.edu"

def fetch_genome_ids():
    # Staphylococcus aureus
    # database query/info
    search_term = '"Staphylococcus aureus"[Organism] AND "complete genome"[Title]'
    db_name = "nucleotide"

    print(f"Querying NCBI with database {db_name} and search term {search_term}...")

    # returns query data as a dictionary
    # retmax = max number of return values
    handle = Entrez.esearch(db=db_name, term=search_term, retmax=127486)
    results = Entrez.read(handle)

    # Count = number of records returned
    count = int(results["Count"])
    print(f"Total genome records found: {count}")
    
    # write all record Ids to file
    file = open("SA_nucleotide_genomes.txt", "w")
    for id in results["IdList"]:
        file.write(str(id) + "\n")
    
    print("Wrote IDs for Staphylococcus aureus to file")

    # Pseudomonas aeruginosa
    # database query/info
    search_term = '"Pseudomonas aeruginosa"[Organism] AND "complete genome"[Title]'
    db_name = "nucleotide"

    print(f"Querying NCBI with database {db_name} and search term {search_term}...")

    # returns query data as a dictionary
    # retmax = max number of return values
    handle = Entrez.esearch(db=db_name, term=search_term, retmax=127486)
    results = Entrez.read(handle)

    # Count = number of records returned
    count = int(results["Count"])
    print(f"Total genome records found: {count}")
    
    # write all record Ids to file
    file = open("PA_nucleotide_genomes.txt", "w")
    for id in results["IdList"]:
        file.write(str(id) + "\n")
    
    print("Wrote IDs for Pseudomonas aeruginosa to file")

fetch_genome_ids()