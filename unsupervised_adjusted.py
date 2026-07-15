import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import pairwise_distances, adjusted_rand_score
import matplotlib.pyplot as plt

SEED = 42
np.random.seed(SEED)

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

def one_hot_positions(series, prefix):
    out = pd.DataFrame(index=series.index)
    for pos in range(3):
        for base in ["A", "C", "G", "T"]:
            col_name = f"{prefix}_pos{pos+1}_{base}"
            out[col_name] = (series.str.upper().str[pos] == base).astype(int)
    return out

def analyze_unsupervised_flexible(input_file, numeric_features, include_codons=True, n_clusters=4, run_label="Run"):
    """
    Runs clustering allowing to pick exactly what traits to use.
    """
    print(f"\n==================== EXECUTION: {run_label} ====================")
    df_full = pd.read_csv(input_file)
    
    # Construct feature matrix based on user choices
    feature_parts = [df_full[numeric_features].reset_index(drop=True)]
    
    if include_codons:
        first3_oh = one_hot_positions(df_full[CSV_COLUMNS["first3"]], "first3")
        last3_oh = one_hot_positions(df_full[CSV_COLUMNS["last3"]], "last3")
        feature_parts.extend([first3_oh.reset_index(drop=True), last3_oh.reset_index(drop=True)])
        
    X = pd.concat(feature_parts, axis=1)
    
    # Scale data (neccessary)
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # K-Means
    kmeans = KMeans(n_clusters=n_clusters, random_state=SEED, n_init=10)
    cluster_labels = kmeans.fit_predict(X_scaled)
    df_full[f'Cluster_{run_label}'] = cluster_labels
    
    # 1. cluster distances!
    # Calculate Euclidean distance between centroids in the scaled feature space
    centroids = kmeans.cluster_centers_
    dist_matrix = pairwise_distances(centroids, metric='euclidean')
    dist_df = pd.DataFrame(dist_matrix, columns=[f"Cluster {i}" for i in range(n_clusters)], index=[f"Cluster {i}" for i in range(n_clusters)])
    
    print("\n[DISTANCE MATRIX] Straight-line distances between cluster centers:")
    print(dist_df.round(2))
    
    #2. decoding the axes (PCA Loadings) 
    pca = PCA(n_components=2, random_state=SEED)
    X_pca = pca.fit_transform(X_scaled)
    
    loadings = pd.DataFrame(
        pca.components_.T, 
        columns=['PC1_Weight', 'PC2_Weight'], 
        index=X.columns
    )
    print("\n[AXIS MEANINGS] Top 5 characteristics shaping the Horizontal Axis (PC1):")
    print(loadings['PC1_Weight'].abs().sort_values(ascending=False).head(5))
    print("\n[AXIS MEANINGS] Top 5 characteristics shaping the Vertical Axis (PC2):")
    print(loadings['PC2_Weight'].abs().sort_values(ascending=False).head(5))
    
    # Plotting
    plt.figure(figsize=(8, 5))
    scatter = plt.scatter(X_pca[:, 0], X_pca[:, 1], c=cluster_labels, cmap='viridis', alpha=0.7, edgecolors='k', linewidths=0.5)
    plt.colorbar(scatter, label='Cluster Group')
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    plt.title(f"Genomic Profile Space: {run_label}")
    plt.grid(True, linestyle=':', alpha=0.5)
    plt.tight_layout()
    plt.show()
    
    return df_full[[CSV_COLUMNS["gene_name"], f'Cluster_{run_label}']]

def compare_cluster_runs(run1_df, run2_df, col_run1, col_run2):
    """
    Compares two different clustering outputs to check alignment.
    """
    # Merge the execution dfs together on the gene name
    merged = pd.merge(run1_df, run2_df, on="gene_name")
    
    print(f"\n==================== COMPARING CONFIGURATIONS ====================")
    # Generate mapping matrix
    matrix = pd.crosstab(merged[col_run1], merged[col_run2])
    print("\n[ALIGNMENT MATRIX] Overlap count between configurations:")
    print(matrix)
    
    # Calculate standard mathematical similarity tracking
    ari_score = adjusted_rand_score(merged[col_run1], merged[col_run2])
    print(f"\nAdjusted Rand Index (ARI) Alignment Score: {ari_score:.4f}")
    if ari_score > 0.75:
        print("--> Excellent stability! The configurations group the genes almost identically.")
    elif ari_score > 0.4:
        print("--> Moderate alignment. Adjusting traits shifted some boundary genes.")
    else:
        print("--> Weak alignment. Changing the inputs completely re-architected the clusters.")


# =====================================================================
# RUN EXPERIMENTS 
# =====================================================================

# EXPERIMENT 1: Cluster using only Regulatory Data (Expression metrics)
features_run1 = [CSV_COLUMNS["log2_fold_change"], CSV_COLUMNS["log_pvalue"]]
run1_results = analyze_unsupervised_flexible(
    input_file="SA_combined.csv", 
    numeric_features=features_run1, 
    include_codons=False, 
    n_clusters=3, 
    run_label="Regulatory_Only"
)

# EXPERIMENT 2: Cluster using Physical/Evolutionary DNA Traits Only
features_run2 = [CSV_COLUMNS["seq_length"], CSV_COLUMNS["gc_pct"]]
run2_results = analyze_unsupervised_flexible(
    input_file="SA_combined.csv", 
    numeric_features=features_run2, 
    include_codons=True, 
    n_clusters=3, 
    run_label="Physical_Sequence_Only"
)

# COMPARE BOTH EXPERIMENTS TO SEE IF THEY MATCH UP
compare_cluster_runs(run1_results, run2_results, "Cluster_Regulatory_Only", "Cluster_Physical_Sequence_Only")