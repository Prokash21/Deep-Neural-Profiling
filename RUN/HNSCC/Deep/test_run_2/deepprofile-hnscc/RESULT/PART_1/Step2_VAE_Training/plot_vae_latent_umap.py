# plot_vae_latent_umap.py
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler
import sys

# Load
file = sys.argv[1]
df = pd.read_csv(file, sep='\t', index_col=0)

# Normalize
scaled = StandardScaler().fit_transform(df)

# t-SNE
tsne = TSNE(n_components=2, random_state=42, perplexity=30)
tsne_result = tsne.fit_transform(scaled)

# Plot
plt.figure(figsize=(8, 6))
sns.scatterplot(x=tsne_result[:, 0], y=tsne_result[:, 1])
plt.title(f"Latent Embedding t-SNE: {file}")
plt.xlabel("t-SNE 1")
plt.ylabel("t-SNE 2")
plt.tight_layout()
plt.savefig(file.replace(".tsv", "_tsne.png"))
plt.savefig(file.replace(".tsv", "_tsne.pdf"))
plt.show()
