import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.manifold import TSNE
import sys

# USAGE: python plot_ensemble_umap.py <embedding_file> [<label_file>]
# Example: python plot_ensemble_umap.py HEAD_NECK_DeepProfile_Training_Embedding_50L.tsv HEAD_NECK_TRAINING_DATA_kmeans_ENSEMBLE_LABELS_50L.txt

embed_file = sys.argv[1]
label_file = sys.argv[2] if len(sys.argv) > 2 else None

df = pd.read_csv(embed_file, sep="\t", index_col=0)

# Use UMAP if available, otherwise t-SNE
try:
    import umap
    reducer = umap.UMAP(random_state=0)
    emb = reducer.fit_transform(df)
    xlab, ylab = "UMAP 1", "UMAP 2"
except ImportError:
    emb = TSNE(n_components=2, random_state=0).fit_transform(df)
    xlab, ylab = "t-SNE 1", "t-SNE 2"

plt.figure(figsize=(8,6))

if label_file:
    labels = pd.read_csv(label_file, sep=",", header=None, names=["cluster"])
    if len(labels) == len(df):
        sns.scatterplot(x=emb[:,0], y=emb[:,1], hue=labels['cluster'].astype(int).astype(str),
                        palette='tab20', s=40, alpha=0.85, legend='full')
        plt.legend(title='Cluster', bbox_to_anchor=(1.01, 1), loc='upper left')
    else:
        print("Label file does not match number of samples—plotting without labels.")
        plt.scatter(emb[:,0], emb[:,1], s=30, alpha=0.8)
else:
    plt.scatter(emb[:,0], emb[:,1], s=30, alpha=0.8)

plt.title("Step 5: Ensemble Embedding UMAP/t-SNE")
plt.xlabel(xlab)
plt.ylabel(ylab)
plt.tight_layout()
plt.savefig("ensemble_embedding_umap.pdf")
print("Saved: ensemble_embedding_umap.pdf")
