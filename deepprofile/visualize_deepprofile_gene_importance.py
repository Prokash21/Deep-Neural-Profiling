#!/usr/bin/env python3
"""
visualize_deepprofile_gene_importance.py

Generates heatmap and barplot for DeepProfile gene importance weights,
and saves PDF and PNG versions without blocking execution.
"""
import os
import pandas as pd
import matplotlib
# Use non-interactive backend so plt.show() won't block
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# ----- CONFIGURATION -----
RESULT_PATH = './ALL_CANCER_FILES/HEAD_NECK/HEAD_NECK_DeepProfile_Ensemble_Gene_Importance_Weights_50L.tsv'
OUTPUT_DIR = './RESULT/PART_1/'
TOP_GENES_HEATMAP = 50  # number of top genes to include in heatmap
TOP_GENES_BAR = 20      # number of top genes to barplot
NODE_FOR_BAR = 0        # latent node index to barplot

# Create output directory if it doesn't exist
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ----- LOAD DATA -----
df = pd.read_csv(RESULT_PATH, sep='\t', index_col=0)

# ----- HEATMAP -----
# select genes by highest variance across latent nodes
gene_variance = df.var(axis=0)
top_genes = gene_variance.sort_values(ascending=False).head(TOP_GENES_HEATMAP).index
heatmap_data = df.loc[:, top_genes]

plt.figure(figsize=(16, 10))
plt.imshow(heatmap_data, aspect='auto', interpolation='nearest')
plt.colorbar(label='Gene Importance', shrink=0.8)
plt.xlabel('Gene', fontsize=14, fontweight='bold')
plt.ylabel('Latent Node', fontsize=14, fontweight='bold')
plt.title('DeepProfile Ensemble Gene Importance (Heatmap)', fontsize=16, fontweight='bold')
plt.xticks(
    ticks=np.arange(len(top_genes)),
    labels=top_genes,
    rotation=90,
    fontsize=8,
    fontname='Arial'
)
plt.yticks(
    ticks=np.arange(df.shape[0]),
    labels=df.index,
    fontsize=8,
    fontname='Arial'
)
plt.tight_layout()
# save to PDF and PNG
heatmap_pdf = os.path.join(OUTPUT_DIR, 'heatmap_final.pdf')
heatmap_png = os.path.join(OUTPUT_DIR, 'heatmap_final.png')
plt.savefig(heatmap_pdf, format='pdf')
plt.savefig(heatmap_png, dpi=300)
plt.close()

# ----- BARPLOT -----
# pick one latent node and its top genes by absolute importance
node_weights = df.iloc[NODE_FOR_BAR].abs().sort_values(ascending=False).head(TOP_GENES_BAR)

plt.figure(figsize=(12, 6))
node_weights.plot(kind='bar')
plt.ylabel('Importance', fontsize=14, fontweight='bold')
plt.xlabel('Gene', fontsize=14, fontweight='bold')
plt.title(f'Top {TOP_GENES_BAR} Genes for Latent Node {NODE_FOR_BAR}', fontsize=16, fontweight='bold')
plt.xticks(rotation=60, fontsize=10, fontname='Arial')
plt.yticks(fontsize=10, fontname='Arial')
plt.tight_layout()
bar_pdf = os.path.join(OUTPUT_DIR, f'top_genes_node_{NODE_FOR_BAR}.pdf')
bar_png = os.path.join(OUTPUT_DIR, f'top_genes_node_{NODE_FOR_BAR}.png')
plt.savefig(bar_pdf, format='pdf')
plt.savefig(bar_png, dpi=300)
plt.close()

print(f"Heatmap saved to: {heatmap_pdf} and {heatmap_png}")
print(f"Barplot saved to: {bar_pdf} and {bar_png}")

# ----- PRINT TOP GENES FOR EACH NODE -----
TOP_GENES_TABLE = 5
print("\nTop genes per latent node:")
for node in df.index:
    top5 = df.loc[node].abs().sort_values(ascending=False).head(TOP_GENES_TABLE).index.tolist()
    print(f"Node {node}: {top5}")
