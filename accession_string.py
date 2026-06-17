import pandas as pd

# takes list of accessions from file and formats it into a submittable query for RGI
def format_string(filename: str):
    df = pd.read_csv(filename)

    res = ""
    i = 0
    max = 25
    for a in df["Accession"]:
        if i >= max:
            break
        
        res += a + ","
        i += 1
    return res

print(format_string("SA_local_ids_summary.csv"))