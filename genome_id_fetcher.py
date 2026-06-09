from Bio import Entrez

Entrez.email = "rcs227@lehigh.edu"

def fetch_genome_ids():
    search_term = '"Staphylococcus aureus"[Organism]'
    db_name = "assembly"

    print(f"Querying NCBI with database {db_name} and search term {search_term}...")

    # find number of genome records
    handle = Entrez.esearch(db=db_name, term=search_term, retmax=127486)
    results = Entrez.read(handle)
    count = int(results["Count"])

    print(f"Total genome records found: {count}")
    
    file = open("genome_ids.txt", "w")
    for id in results["IdList"]:
        file.write(str(id) + "\n")
    
    print("Wrote IDs to file")

fetch_genome_ids()