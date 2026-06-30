import sys
import pandas as pd

if len(sys.argv) < 3:
    print("Error: Expected 2-3 filenames")
    print("Usage: python3 file1 file2 output_file")
    sys.exit(1)

# combines the csv file with regulation data with the csv file with at/gc composition data we found
def combine(reg_file, comp_file, output_file="output.csv"):
    reg_df = pd.read_csv(reg_file)
    comp_df = pd.read_csv(comp_file)

    # columns in composition file to exclude from output file
    exclusions = ["Genome_ID", "Genome_Name"]

    comp_df_filtered = comp_df.drop(columns=exclusions)

    combined_df = pd.merge(reg_df, comp_df_filtered, left_on=["gene_name"], right_on=["Gene_Found"], how='inner')
    combined_df = combined_df.drop(columns="Gene_Found")
    combined_df.to_csv(output_file, index = False)


if len(sys.argv) >= 4:
    combine(sys.argv[1], sys.argv[2], sys.argv[3])
else:
    combine(sys.argv[1], sys.argv[2])