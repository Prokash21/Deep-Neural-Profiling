import pandas as pd
import matplotlib.pyplot as plt
import sys
from matplotlib import rcParams
import os
# USAGE: python plot_top_genes_per_latent.py <attribution_file> <latent_node_number>
# Example: python plot_top_genes_per_latent.py HEAD_NECK_DeepProfile_Ensemble_Weights_50L.tsv 0


# === PDF + FONT SETTINGS ===
rcParams['pdf.fonttype'] = 42       # Illustrator-compatible text
rcParams['font.family'] = 'Arial'   # Use clean, editable Arial font


att_file = sys.argv[1]
node = int(sys.argv[2])

df = pd.read_csv(att_file, sep="\t", index_col=0)
# Get top 10 positive and top 10 negative genes (absolute attribution)
vals = df.iloc[node]
top_genes = vals.abs().sort_values(ascending=False).head(10).index
top_vals = vals.loc[top_genes]

plt.figure(figsize=(7,4))
top_vals.plot(kind="bar", rot=4)
plt.title(f"Top 10 Genes for Latent Node {node}")
plt.ylabel("Attribution Score")
plt.tight_layout()
plt.savefig(f"top_genes_latent_{node}.pdf")
print(f"Saved: top_genes_latent_{node}.pdf")
