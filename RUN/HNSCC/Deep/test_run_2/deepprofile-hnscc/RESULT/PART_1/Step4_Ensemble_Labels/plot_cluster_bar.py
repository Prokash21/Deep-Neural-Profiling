import pandas as pd
import matplotlib.pyplot as plt
import sys

# USAGE: python plot_cluster_bar.py <label_file>
# Example: python plot_cluster_bar.py HEAD_NECK_TRAINING_DATA_kmeans_ENSEMBLE_LABELS_50L.txt

label_file = sys.argv[1]
labels = pd.read_csv(label_file, sep=",", header=None, names=["cluster"])

sizes = labels['cluster'].astype(int).value_counts().sort_index()

plt.figure(figsize=(7,4))
plt.bar(sizes.index.astype(str), sizes.values)
plt.xlabel("Cluster")
plt.ylabel("# of samples")
plt.title("Cluster Size Distribution")
plt.tight_layout()
plt.savefig("ensemble_cluster_sizes.pdf")
print("Saved: ensemble_cluster_sizes.pdf")
