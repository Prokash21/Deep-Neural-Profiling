import pandas as pd
import matplotlib.pyplot as plt
import sys

# Usage: python plot_ig_bar_top.py <IG_file.tsv> <latent_idx> <fold>
# Example: python plot_ig_bar_top.py HEAD_NECK_DATA_VAE_Cluster_Weights_TRAINING_50L_fold0.tsv 0 0

if len(sys.argv) < 4:
    print("Usage: python plot_ig_bar_top.py <IG_file.tsv> <latent_idx> <fold>")
    sys.exit(1)

file = sys.argv[1]
latent_idx = int(sys.argv[2])
fold = sys.argv[3]

df = pd.read_csv(file, sep="\t", index_col=0)
top_n = 10
abs_scores = df.iloc[latent_idx].abs().sort_values(ascending=False)
top = abs_scores[:top_n]

plt.figure(figsize=(8, 3))
plt.bar(top.index.astype(str), top.values)
plt.xlabel("PCA Component")
plt.ylabel("IG Score")
plt.title(f"Top PCA for Latent {latent_idx} (Fold {fold})")
plt.tight_layout()
plt.savefig(f"IG_bar_latent{latent_idx}_fold{fold}.pdf")
