# NOTE: Unsupervised learning maps out the natural hidden architecture of the genes to 
# group them by overall profile similarities (so like expression patterns + physical
# sequence metrics combined). As such, the approach below uses K-Means Clustering (to group
# the genes and PCA (Principal Component Analysis) (to compress the features down to 2D space
# in order to visibly see the clusters). 

# This model is exploratory -- it determines whether the organism has distinct sub-classes
# or evolutionary profiles of genes. Also, since there arent any hidden labels, there doesnt
# need to be separate validation slices -- every metric acts as a cooperative feature for 
# structural grouping. 

# NOTE: the resulting output file (SA_unsupervised_clustered.csv) will append a literal cluster
# signature tag directly onto each gene row. This can then be used to filter the file by cluster
# number to locate specific functional categories or operons that naturally adjust or assemble in uniform ways.

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

# Set random seed for reproducibility
SEED = 42
np.random.seed(SEED)

# Universal column dictionary (same as NN model for continuity)
CSV_COLUMNS = {
    "gene_name": "gene_name",
    "log2_fold_change": "log2FoldChange",
    "pvalue": "pvalue",
    "log_pvalue": "log_pvalue",
    "seq_length": "Sequence_Length",
    "first3": "First_3",
    "last3": "Last_3",
    "gc_pct": "GC_pct",
    "at_pct": "AT_pct",
}

BASES = ["A", "C", "G", "T"]

def one_hot_positions(series, prefix):
    """One-hot encodes 3-letter codon structures into 12 distinct 
    binary markers, matching your colleague's exact feature pipeline."""
    out = pd.DataFrame(index=series.index)
    for pos in range(3):
        for base in BASES:
            col_name = f"{prefix}_pos{pos+1}_{base}"
            out[col_name] = (series.str.upper().str[pos] == base).astype(int)
    return out


def run_unsupervised_analysis(input_file, n_clusters=4):
    print(f"--- Loading Dataset: {input_file} ---")
    df_full = pd.read_csv(input_file)
    
    # 1. Process categorical sequence columns using the baseline pipeline
    first3_oh = one_hot_positions(df_full[CSV_COLUMNS["first3"]], "first3")
    last3_oh = one_hot_positions(df_full[CSV_COLUMNS["last3"]], "last3")
    
    # 2. Select numerical features to cluster on
    # In unsupervised learning, nothing is held back as a hidden "target"! 
    # so, include both regulatory variables and physiological sequence features.
    numeric_cols = [
        CSV_COLUMNS["log2_fold_change"],
        CSV_COLUMNS["log_pvalue"],
        CSV_COLUMNS["seq_length"],
        CSV_COLUMNS["gc_pct"]
    ]
    
    # 3. Concatenate all variables into a single feature array
    X = pd.concat([
        df_full[numeric_cols].reset_index(drop=True), 
        first3_oh.reset_index(drop=True), 
        last3_oh.reset_index(drop=True)
    ], axis=1)
    
    print(f"Generated Feature Profile Grid with shape: {X.shape}")
    
    # 4. Standard Scaling (MANDATORY for K-Means and PCA)
    # Because Sequence Lengths are in the hundreds/thousands while GC fractions 
    # and Log2FC are tiny, unscaled profiles would cluster entirely on sequence length.
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # 5. Execute K-Means Clustering Algorithm
    print(f"Running K-Means Clustering (Targeting K={n_clusters} modules)...")
    kmeans = KMeans(n_clusters=n_clusters, random_state=SEED, n_init=10)
    cluster_labels = kmeans.fit_predict(X_scaled)
    
    # Assign the resulting signatures back to the original layout
    df_full['Cluster_Assignment'] = cluster_labels
    
    # 6. Execute Dimensionality Reduction (PCA) for Graphical Plotting
    print("Executing Principal Component Analysis for 2D Projection...")
    pca = PCA(n_components=2, random_state=SEED)
    X_pca = pca.fit_transform(X_scaled)
    
    # Read variance explanations to document mathematically
    var_ratios = pca.explained_variance_ratio_
    print(f"  - PC1 Explains: {var_ratios[0]*100:.2f}% of data variance")
    print(f"  - PC2 Explains: {var_ratios[1]*100:.2f}% of data variance")
    
    # 7. Generate Cluster Scattering Map
    plt.figure(figsize=(10, 6))
    scatter = plt.scatter(
        X_pca[:, 0], X_pca[:, 1], 
        c=cluster_labels, 
        cmap='viridis', 
        alpha=0.8, 
        edgecolors='black', 
        linewidths=0.5
    )
    
    plt.colorbar(scatter, label='Assigned Cluster Group')
    plt.xlabel(f"Principal Component 1 ({var_ratios[0]*100:.1f}%)")
    plt.ylabel(f"Principal Component 2 ({var_ratios[1]*100:.1f}%)")
    plt.title(f"Unsupervised Genomic Map: Staphylococcus aureus Genes (Clusters: {n_clusters})")
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.tight_layout()
    plt.show()
    
    # Print quantitative breakdown profiles for documentation
    print("\n=== Cluster Gene Allocation Count ===")
    print(df_full['Cluster_Assignment'].value_counts().sort_index())
    
    # Save the clustered data mapping out to a fresh file
    output_filename = "SA_unsupervised_clustered.csv"
    df_full.to_csv(output_filename, index=False)
    print(f"\nSaved master cluster configuration mapping to '{output_filename}'")
    
    return df_full

# --- Execution ---
# Call the function directly on the target combined table file
clustered_df = run_unsupervised_analysis("SA_combined.csv", n_clusters=4)