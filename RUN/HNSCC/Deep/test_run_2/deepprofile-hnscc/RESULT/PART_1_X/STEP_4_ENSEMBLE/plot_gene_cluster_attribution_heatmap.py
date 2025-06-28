import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

labels_file = "HEAD_NECK_TRAINING_DATA_kmeans_ENSEMBLE_LABELS_50L.txt"
ig_file0 = "HEAD_NECK_DATA_VAE_Cluster_Weights_TRAINING_50L_fold0.tsv"
ig_file1 = "HEAD_NECK_DATA_VAE_Cluster_Weights_TRAINING_50L_fold1.tsv"
out_file = "gene_cluster_attribution_heatmap.pdf"

# 1. Average the two IG files
ig0 = pd.read_csv(ig_file0, sep='\t', index_col=0)
ig1 = pd.read_csv(ig_file1, sep='\t', index_col=0)
ig = (ig0 + ig1) / 2.0

# 2. Read cluster assignments (single column, no header)
labels = pd.read_csv(labels_file, sep="\t", header=None, names=["cluster"])

# 3. Make sure both have the same number of latents/columns
n_latents = min(len(labels), ig.shape[1])
if len(labels) != ig.shape[1]:
    print(f"WARNING: label file has {len(labels)} rows, IG matrix has {ig.shape[1]} columns. Using first {n_latents} of each.")
labels = labels.iloc[:n_latents]
ig = ig.iloc[:, :n_latents]

# 4. Get top N genes
N = 20
gene_scores = ig.abs().sum(axis=1)
top_genes = gene_scores.nlargest(N).index
ig_top = ig.loc[top_genes]

# 5. Compute mean IG per cluster using only **valid** latent indices
cluster_ids = np.sort(labels['cluster'].unique())
cluster_ig = np.zeros((len(top_genes), len(cluster_ids)))
for j, cluster_id in enumerate(cluster_ids):
    # These are latent indices in the label file (must be less than n_latents)
    latent_indices = labels.index[labels['cluster'] == cluster_id].tolist()
    # Filter for indices that are in-bounds for IG matrix
    latent_indices = [idx for idx in latent_indices if idx < ig_top.shape[1]]
    if len(latent_indices) > 0:
        cluster_ig[:, j] = ig_top.iloc[:, latent_indices].mean(axis=1).values
    else:
        # If no valid latents for this cluster, fill with NaN
        cluster_ig[:, j] = np.nan

# 6. Plot
plt.figure(figsize=(12, 6))
im = plt.imshow(cluster_ig, aspect='auto', origin='lower', cmap='coolwarm')
plt.colorbar(im, label="Mean IG Attribution")
plt.yticks(np.arange(len(top_genes)), top_genes)
plt.xticks(np.arange(len(cluster_ids)), cluster_ids.astype(int), rotation=90)
plt.xlabel("Ensemble Cluster ID")
plt.ylabel("Gene")
plt.title(f"Gene Attribution Patterns Across Ensemble Clusters (Top {N} Genes)")
plt.tight_layout()
plt.savefig(out_file)
plt.close()
print("Saved gene-cluster attribution heatmap:", out_file)
