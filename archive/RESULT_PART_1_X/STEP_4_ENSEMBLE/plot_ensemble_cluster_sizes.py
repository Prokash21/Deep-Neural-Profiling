import pandas as pd
import matplotlib.pyplot as plt

labels_file = "HEAD_NECK_TRAINING_DATA_kmeans_ENSEMBLE_LABELS_50L.txt"
out_file = "ensemble_cluster_sizes_50L.pdf"

# Read the single-column cluster assignments
clusters = pd.read_csv(labels_file, sep="\t", header=None)[0]

# Latent IDs are just their index (0, 1, 2, ...)
latents = pd.Series(range(len(clusters)), name="latent")
labels = pd.DataFrame({"latent": latents, "cluster": clusters})

cluster_counts = labels['cluster'].value_counts().sort_index()

plt.figure(figsize=(9, 4))
plt.bar(cluster_counts.index.astype(int), cluster_counts.values, color='#4477AA')
plt.xlabel("Ensemble Cluster ID")
plt.ylabel("Number of Latents")
plt.title("Cluster Sizes for Ensemble Labels (K=50)")
plt.tight_layout()
plt.savefig(out_file)
plt.close()
print("Saved bar plot:", out_file)
