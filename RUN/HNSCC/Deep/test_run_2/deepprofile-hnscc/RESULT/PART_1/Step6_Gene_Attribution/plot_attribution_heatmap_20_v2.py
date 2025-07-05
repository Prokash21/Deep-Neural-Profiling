import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib.gridspec as gridspec

# === PDF + FONT SETTINGS ===
rcParams['pdf.fonttype'] = 42       # Illustrator-compatible text
rcParams['font.family'] = 'Arial'   # Use clean, editable Arial font

# ── CONFIG ──────────────────────────────────────────────
ig_file = "HEAD_NECK_DeepProfile_Ensemble_Weights_50L.tsv"
out_file = "ensemble_ig_top20_dendro_heatmap.pdf"  # Output PDF saved in working directory

# ── LOAD FINAL ENSEMBLED ATTRIBUTION MATRIX ─────────────
ig = pd.read_csv(ig_file, sep="\t", index_col=0)  # shape: latent × gene

# TRANSPOSE: rows = genes, columns = latents
ig = ig.T

# ── PICK TOP 20 GENES ───────────────────────────────────
gene_scores = ig.abs().sum(axis=1)
top20 = gene_scores.nlargest(20).index
mat = ig.loc[top20].T

# ── HIERARCHICAL CLUSTERING ─────────────────────────────
row_link = linkage(mat.values, method="average", metric="correlation")
col_link = linkage(mat.values.T, method="average", metric="correlation")

row_dendr = dendrogram(row_link, no_plot=True)
col_dendr = dendrogram(col_link, no_plot=True)
row_idx, col_idx = row_dendr["leaves"], col_dendr["leaves"]
mat_ord = mat.iloc[row_idx, :].iloc[:, col_idx]

# ── DRAW FIGURE ─────────────────────────────────────────
fig = plt.figure(figsize=(5, 8))
gs = gridspec.GridSpec(2, 2,
    width_ratios=(0.15, 4),
    height_ratios=(1, 4),
    wspace=0.0, hspace=0.0)

# Column dendrogram (top)
ax_col = fig.add_subplot(gs[0,1])
dendrogram(col_link, ax=ax_col, color_threshold=None,
           labels=mat_ord.columns, orientation='top')
ax_col.axis('off')

# Row dendrogram (left)
ax_row = fig.add_subplot(gs[1,0])
dendrogram(row_link, ax=ax_row, color_threshold=None,
           labels=mat_ord.index, orientation='left')
ax_row.axis('off')

# Heatmap
ax_heat = fig.add_subplot(gs[1,1])
im = ax_heat.matshow(mat_ord.values, aspect='auto', origin='lower', cmap='viridis')

# Axis labels
ax_heat.set_xticks(np.arange(len(mat_ord.columns)))
ax_heat.set_xticklabels(mat_ord.columns, rotation=90, fontsize=9)
ax_heat.xaxis.set_ticks_position('bottom')
ax_heat.tick_params(axis='x', labelbottom=True, labeltop=False)

ax_heat.set_yticks(np.arange(len(mat_ord.index)))
ax_heat.set_yticklabels(mat_ord.index, fontsize=9)
ax_heat.yaxis.set_ticks_position('right')
ax_heat.tick_params(axis='y', labelleft=False, labelright=True)

# Colorbar
cax = fig.add_axes([0.965, 0.15, 0.015, 0.7])
fig.colorbar(im, cax=cax, label='Mean Attribution Score')

# Title and save
plt.suptitle("Hierarchically Clustered Heatmap\nTop 20 Genes × 50 Latents (Ensembled)", y=0.98)
plt.tight_layout(rect=[0,0,0.96,1])

plt.savefig(out_file)
print("Saved dendrogram heatmap to", out_file)
