import pandas as pd

def assign_regulation(val):
    # Handle any potential blank/NaN values from the map
    if pd.isna(val):
        return 0 
        
    if val >= 1:
        return 1
    elif val <= -1:
        return -1
    else:
        return 0

def combine(shine_file, dna_file, output_file):
    shine_df = pd.read_csv(shine_file)
    dna_df = pd.read_csv(dna_file)

    unique_dna = dna_df.drop_duplicates(subset="gene_name")


    shine_df['log2FoldChange'] = shine_df['Gene_Name'].map(unique_dna.set_index("gene_name")['log2FoldChange'])

    shine_df['regulation'] = shine_df['log2FoldChange'].apply(assign_regulation)
    
    # Save the result
    shine_df.to_csv(output_file, index=False)
    print(f"Successfully saved to {output_file}")

combine("SA_shine.csv", "SA_DNA.csv", "SA_shine_data.csv")