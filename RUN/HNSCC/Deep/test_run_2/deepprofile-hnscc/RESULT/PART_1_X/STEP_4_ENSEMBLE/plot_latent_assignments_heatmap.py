import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

labels_file = "HEAD_NECK_TRAINING_DATA_kmeans_ENSEMBLE_LABELS_50L.txt"
out_file = "latent_assignments_heatmap_50L.pdf"

clusters = pd.read_csv(labels_file, sep="\t", header=None)[0]
latents = pd.Series(range(len(clusters)), name="latent")
labels = pd.DataFrame({"latent": latents, "cluster": clusters})

num_latents = len(labels)
cluster_ids = np.sort(labels['cluster'].unique())
num_clusters = len(cluster_ids)

assignment = np.zeros((num_latents, num_clusters))
for i, cluster_id in enumerate(labels['cluster']):
    j = np.where(cluster_ids == cluster_id)[0][0]
    assignment[i, j] = 1

plt.figure(figsize=(10, 8))
plt.imshow(assignment, aspect='auto', cmap='Blues', origin='lower')
plt.title("Latent-to-Ensemble Cluster Assignment (50 Latents × K Clusters)")
plt.xlabel("Ensemble Cluster ID")
plt.ylabel("Latent Feature ID")
plt.xticks(np.arange(num_clusters), cluster_ids.astype(int), rotation=90)
plt.yticks(np.arange(num_latents), labels['latent'])
plt.tight_layout()
plt.savefig(out_file)
plt.close()
print("Saved assignment heatmap:", out_file)
