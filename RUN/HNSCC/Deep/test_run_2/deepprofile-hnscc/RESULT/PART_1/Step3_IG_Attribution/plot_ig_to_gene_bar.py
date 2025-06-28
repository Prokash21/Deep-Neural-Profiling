import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

# Usage: python plot_ig_to_gene_bar.py <IG_file.tsv> <latent_idx> <fold> <COMPONENTS_file.tsv>
# Example: python plot_ig_to_gene_bar.py HEAD_NECK_DATA_VAE_Cluster_Weights_TRAINING_50L_fold0.tsv 0 0 HEAD_NECK_DATA_TOP2_JOINED_PCA_500L_COMPONENTS.tsv

if len(sys.argv) < 5:
    print("Usage: python plot_ig_to_gene_bar.py <IG_file.tsv> <latent_idx> <fold> <COMPONENTS_file.tsv>")
    sys.exit(1)

ig_file = sys.argv[1]
latent_idx = int(sys.argv[2])
fold = sys.argv[3]
comp_file = sys.argv[4]

# IG file
ig = pd.read_csv(ig_file, sep="\t", index_col=0)
# PCA components
comp = pd.read_csv(comp_file, sep="\t", index_col=0)

# Get top PC components for latent
top_pcs = ig.iloc[latent_idx].abs().sort_values(ascending=False)[:5].index.astype(int)
# For each PC, get top genes
top_gene_set = set()
for pc in top_pcs:
    top_gene_set.update(comp.iloc[:, pc].abs().sort_values(ascending=False)[:10].index)
# Build heatmap data
heat = comp.loc[list(top_gene_set), [str(pc) for pc in top_pcs]]
plt.figure(figsize=(8, 4))
sns.heatmap(heat, cmap="coolwarm", center=0)
plt.title(f"Genes for Latent {latent_idx} via Top PCs (Fold {fold})")
plt.tight_layout()
plt.savefig(f"IG_latent{latent_idx}_to_gene_fold{fold}.pdf")
