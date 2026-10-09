import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys

# USAGE: python plot_ensemble_heatmap.py <embedding_file>
# Example: python plot_ensemble_heatmap.py HEAD_NECK_DeepProfile_Training_Embedding_50L.tsv

embed_file = sys.argv[1]
df = pd.read_csv(embed_file, sep="\t", index_col=0)

plt.figure(figsize=(14,8))
sns.heatmap(df, cmap="viridis", yticklabels=False, cbar_kws={'label': 'Latent Intensity'})
plt.title("Step 5: Ensemble Embedding Heatmap")
plt.xlabel("Latent Feature")
plt.ylabel("Sample")
plt.tight_layout()
plt.savefig("ensemble_embedding_heatmap.pdf")
print("Saved: ensemble_embedding_heatmap.pdf")
