import pandas as pd
import matplotlib.pyplot as plt
import os
from matplotlib import rcParams

# === PDF + FONT SETTINGS ===
rcParams['pdf.fonttype'] = 42       # Illustrator-compatible text
rcParams['font.family'] = 'Arial'   # Use clean, editable Arial font

# CONFIG
att_file = "HEAD_NECK_DeepProfile_Ensemble_Weights_50L.tsv"  # ensemble attribution file
output_file = "top_genes_all_latents.pdf"
n_latents = 50

# Load attribution data
df = pd.read_csv(att_file, sep="\t", index_col=0)  # shape: latent × gene

# Create figure with subplots
fig, axes = plt.subplots(nrows=10, ncols=5, figsize=(18, 24))
axes = axes.flatten()

# Iterate through each latent node
for node in range(n_latents):
    ax = axes[node]
    vals = df.iloc[node]
    top_genes = vals.abs().sort_values(ascending=False).head(10).index
    top_vals = vals.loc[top_genes]
    
    top_vals.plot(kind="bar", ax=ax)
    ax.set_title(f"Latent {node}", fontsize=9)
    ax.tick_params(axis='x', labelrotation=90, labelsize=7)
    ax.tick_params(axis='y', labelsize=8)
    ax.set_ylabel("Score", fontsize=7)

# Adjust layout
plt.suptitle("Top 10 Genes for Each Latent Node", fontsize=16)
plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.savefig(output_file)
print(f"Saved: {output_file}")
