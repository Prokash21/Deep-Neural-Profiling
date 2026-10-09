import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys

# USAGE: python plot_ensemble_cluster_means_heatmap.py <embedding_file> <label_file>
embed_file = sys.argv[1]
label_file = sys.argv[2]

df = pd.read_csv(embed_file, sep="\t", index_col=0)
labels = pd.read_csv(label_file, sep=",", header=None, names=["cluster"])
labels.index = df.index
df['cluster'] = labels['cluster'].astype(int)

cluster_means = df.groupby('cluster').mean().sort_index()

plt.figure(figsize=(14,7))
sns.heatmap(cluster_means, cmap="viridis", cbar_kws={'label': 'Mean Latent Intensity'})
plt.title("Step 5: Cluster-Averaged Ensemble Embedding Heatmap")
plt.xlabel("Latent Feature")
plt.ylabel("Cluster")
plt.tight_layout()
plt.savefig("ensemble_embedding_cluster_means_heatmap.pdf")
print("Saved: ensemble_embedding_cluster_means_heatmap.pdf")
