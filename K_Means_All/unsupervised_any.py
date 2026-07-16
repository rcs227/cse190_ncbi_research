import os
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import pairwise_distances, adjusted_rand_score
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_score

SEED = 42
np.random.seed(SEED)

def one_hot_positions(series, prefix):
    """
    Encodes 3-letter codon structures (like ATG, TAA) into binary markers.
    """
    out = pd.DataFrame(index=series.index)
    for pos in range(3):
        for base in ["A", "C", "G", "T"]:
            col_name = f"{prefix}_pos{pos+1}_{base}"
            out[col_name] = (series.astype(str).str.upper().str[pos] == base).astype(int)
    return out


def load_and_prepare_data(input_files, default_id="ID"):
    """
    Flexible loader that handles:
      1. A single file path (string).
      2. Multiple files with the exact same ID column name.
      3. Multiple files with different ID column names, merging them together.
    
    Format for input_files:
      - Single file: "file1.csv"
      - Multi-file (same ID name): ["file1.csv", "file2.csv"]
      - Multi-file (different ID names): [("file1.csv", "gene_name"), ("file2.csv", "Target_Searched")]
    
    NOTE: Stacks files vertically (concatenation) so you cluster ALL genes 
    from both files together.
    
    """
    if isinstance(input_files, str):
        print(f"Loading single dataset: {input_files}")
        return pd.read_csv(input_files), default_id

    elif isinstance(input_files, list):
        print(f"Concatenating (stacking) {len(input_files)} datasets vertically...")
        dfs = []
        for item in input_files:
            # Handle both string paths and (path, id_col) tuples
            if isinstance(item, tuple):
                path, id_col = item
                df = pd.read_csv(path).rename(columns={id_col: default_id})
            else:
                path = item
                df = pd.read_csv(path)
                
            dfs.append(df)
            
        # Stack vertically, aligning columns with the same headers
        master_df = pd.concat(dfs, axis=0, ignore_index=True)
        return master_df, default_id            
    else:
        raise TypeError("input_files must be a string path or a list of files/tuples.")


def find_optimal_clusters(X_scaled, max_k=10, run_label="Run"):
    """
    Computes Inertia and Silhouette scores for K = 2 to max_k.
    Plots the results and automatically returns the best K.
    """
    inertias = []
    silhouette_scores = []
    k_range = range(2, max_k + 1)
    
    print(f"\nEvaluating optimal cluster count for '{run_label}' (Testing K = 2 to {max_k})...")
    
    for k in k_range:
        kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
        labels = kmeans.fit_predict(X_scaled)
        
        # Calculate Elbow Metric (Inertia)
        inertias.append(kmeans.inertia_)
        
        # Calculate Silhouette Metric
        sil_avg = silhouette_score(X_scaled, labels)
        silhouette_scores.append(sil_avg)
        print(f"  -> K = {k} | Avg Silhouette Score: {sil_avg:.4f} | Inertia: {inertias[-1]:.2f}")
        
    # AUTOMATIC SELECTION: Find K with the highest Silhouette Score
    best_k_index = np.argmax(silhouette_scores)
    best_k = k_range[best_k_index]
    print(f"\n[AUTO-SELECTION] Best K mathematically is **K = {best_k}** with Silhouette Score: {silhouette_scores[best_k_index]:.4f}")
    
    # Optional: Plot the curves and save them as an image so you can inspect them
    fig, ax1 = plt.subplots(figsize=(10, 5))
    
    # Plot Inertia (Elbow)
    color = 'tab:blue'
    ax1.set_xlabel('Number of Clusters (K)')
    ax1.set_ylabel('Inertia (Lower is better)', color=color)
    ax1.plot(k_range, inertias, marker='o', color=color, linestyle='-', label='Inertia')
    ax1.tick_params(axis='y', labelcolor=color)
    
    # Plot Silhouette Scores
    ax2 = ax1.twinx()  
    color = 'tab:red'
    ax2.set_ylabel('Silhouette Score (Higher is better)', color=color)
    ax2.plot(k_range, silhouette_scores, marker='s', color=color, linestyle='--', label='Silhouette')
    ax2.tick_params(axis='y', labelcolor=color)
    
    plt.title(f'Finding Optimal K: {run_label}')
    fig.tight_layout()
    
    image_filename = f"optimization_plot_{run_label}.png"
    plt.savefig(image_filename, dpi=300)
    plt.close()
    print(f"[PLOT] Saved diagnostic curves to '{image_filename}'")
    
    return best_k


def analyze_unsupervised_universal(input_files, id_column, numeric_features, 
                                   codon_features=None, categorical_features=None, 
                                   n_clusters=4, run_label="Run"):
    """
    A completely universal clustering tool for any CSV structure.
    """
    print(f"\n==================== EXECUTION: {run_label} ====================")
    
    # 1. Load and dynamically merge files
    df_full, unified_id_col = load_and_prepare_data(input_files, default_id=id_column)
    print(f"Dataset successfully constructed with shape: {df_full.shape}")
    
    # 2. Extract numeric columns
    feature_parts = [df_full[numeric_features].reset_index(drop=True)]
    
    # 3. Handle 3-letter biological codons if requested
    if codon_features:
        for codon_col in codon_features:
            codon_oh = one_hot_positions(df_full[codon_col], codon_col.lower())
            feature_parts.append(codon_oh.reset_index(drop=True))
            
    # 4. Handle any other generic categorical columns using standard one-hot encoding
    if categorical_features:
        for cat_col in categorical_features:
            cat_oh = pd.get_dummies(df_full[cat_col], prefix=cat_col.lower(), dtype=int)
            feature_parts.append(cat_oh.reset_index(drop=True))
            
    # Combine everything into the feature matrix X
    X = pd.concat(feature_parts, axis=1)
    
    # Clean any missing/NaN values by filling them with the column mean
    X = X.fillna(X.mean())
    
    # 5. Scale the features (mandatory)
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # 6. If no cluster number is provided, let the computer find the best one!
    if n_clusters is None:
        n_clusters = find_optimal_clusters(X_scaled, max_k=10, run_label=run_label)

    # 7. Now run K-Means using the automated cluster count
    print(f"Running final K-Means Clustering using selected K = {n_clusters}...")
    kmeans = KMeans(n_clusters=n_clusters, random_state=SEED, n_init=10)
    cluster_labels = kmeans.fit_predict(X_scaled)
    
    # Add cluster assignments back to the master dataframe
    cluster_col_name = f'Cluster_{run_label}'
    df_full[cluster_col_name] = cluster_labels
    
    # 7. EXPORT INDIVIDUAL CSV FILES FOR EACH CLUSTER
    print("\n[EXPORT] Splitting dataset into separate cluster files:")
    for cluster_id in range(n_clusters):
        cluster_data = df_full[df_full[cluster_col_name] == cluster_id]
        output_path = f"{run_label}_cluster_{cluster_id}.csv"
        cluster_data.to_csv(output_path, index=False)
        print(f"  --> Saved {len(cluster_data)} rows to '{output_path}'")
    
    # 8. Calculate Cluster Distances
    centroids = kmeans.cluster_centers_
    dist_matrix = pairwise_distances(centroids, metric='euclidean')
    dist_df = pd.DataFrame(
        dist_matrix, 
        columns=[f"Cluster {i}" for i in range(n_clusters)], 
        index=[f"Cluster {i}" for i in range(n_clusters)]
    )
    print("\n[DISTANCE MATRIX] Straight-line distances between cluster centers:")
    print(dist_df.round(2))
    
    # 9. Decode PCA Axis Meanings
    pca = PCA(n_components=2, random_state=SEED)
    X_pca = pca.fit_transform(X_scaled)
    
    loadings = pd.DataFrame(
        pca.components_.T, 
        columns=['PC1_Weight', 'PC2_Weight'], 
        index=X.columns
    )
    print("\n[AXIS MEANINGS] Top 3 characteristics shaping Horizontal Axis (PC1):")
    print(loadings['PC1_Weight'].abs().sort_values(ascending=False).head(3))
    print("\n[AXIS MEANINGS] Top 3 characteristics shaping Vertical Axis (PC2):")
    print(loadings['PC2_Weight'].abs().sort_values(ascending=False).head(3))
    
    # 10. Generate Plot and Save as an Image (To avoid terminal hanging)
    plt.figure(figsize=(8, 5))
    scatter = plt.scatter(X_pca[:, 0], X_pca[:, 1], c=cluster_labels, cmap='viridis', alpha=0.7, edgecolors='k', linewidths=0.5)
    plt.colorbar(scatter, label='Cluster Group')
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    plt.title(f"Genomic Profile Space: {run_label}")
    plt.grid(True, linestyle=':', alpha=0.5)
    plt.tight_layout()
    
    # Save image
    image_filename = f"plot_{run_label}.png"
    plt.savefig(image_filename, dpi=300)
    plt.close() # Close memory immediately
    print(f"\n[PLOT] Saved visualization to '{image_filename}'.")
    
    return df_full[[unified_id_col, cluster_col_name]]

# =====================================================================
# RUN EXPERIMENTS
# =====================================================================

# --- EXAMPLE 1: Run on your already combined dataset ---
# run1_results = analyze_unsupervised_universal(
#     input_files="SA_combined.csv", 
#     id_column="gene_name",  # The column to keep as identifier
#     numeric_features=["log2FoldChange", "log_pvalue", "Sequence_Length", "GC_pct"],
#     codon_features=["First_3", "Last_3"], # Uses the biological codon encoding
#     n_clusters=4, 
#     run_label="Single_File_Combined_Run"
# )


# --- EXAMPLE 2: Merge raw files directly (Automatic Combine + Cluster!) ---
# This bypasses needing to run 'combine_csv.py' beforehand!
# It loads both, renames 'gene_name' and 'Target_Searched' to a matching ID name, merges them, and clusters.
raw_files_with_different_headers = [
    ("K_Means_All/CDMAdevsCDMRib.csv", "gene_name"), 
    ("K_Means_All/M9AdevsM9Rib.csv", "gene_name")
]

run2_results = analyze_unsupervised_universal(
    input_files=raw_files_with_different_headers, 
    id_column="Unified_Gene_ID", # Merges on this standard name
    #numeric_features=["log2FoldChange", "log_pvalue", "Sequence_Length", "GC_pct", "AT_pct"],
    numeric_features=["log2FoldChange", "pvalue", "log pvalue"],
    #codon_features=["First_3", "Last_3"],
    #n_clusters=3, <-- ONLY INCLUDE IF YOU WANT TO HARDCODE NUM OF CLUSTERS
    run_label="Direct_Merged_Run"
)


# --- EXAMPLE 3: Run on completely different CSV files (e.g. non-biological or general data) ---
# If you have general data (like car stats, user metrics, etc.), you can do:
# run3_results = analyze_unsupervised_universal(
#     input_files="car_dataset.csv",
#     id_column="car_model_name",
#     numeric_features=["horsepower", "miles_per_gallon", "weight"],
#     categorical_features=["fuel_type", "transmission_type"], # Generic one-hotting!
#     n_clusters=5,
#     run_label="Car_Clustering"
# )