import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram

# ── CONFIG ──────────────────────────────────────────────
base    = r"E:\DWCT\Papia Mam\Cancer\GitHub\Cancer_Project_BMU\RUN\HNSCC\Deep\test_run_2\deepprofile-hnscc"
ig_dir  = os.path.join(base, "ALL_CANCER_FILES", "HEAD_NECK", "VAE_FILES")
out_dir = os.path.join(base, "RESULT", "PART_1", "STEP_3_IG")
os.makedirs(out_dir, exist_ok=True)

# ── LOAD & AVERAGE ──────────────────────────────────────
f0 = os.path.join(ig_dir, "HEAD_NECK_DATA_VAE_Cluster_Weights_TRAINING_50L_fold0.tsv")
f1 = os.path.join(ig_dir, "HEAD_NECK_DATA_VAE_Cluster_Weights_TRAINING_50L_fold1.tsv")
df0 = pd.read_csv(f0, sep="\t", index_col=0)
df1 = pd.read_csv(f1, sep="\t", index_col=0)
ig  = (df0 + df1) / 2.0

# ── PICK TOP 20 GENES ────────────────────────────────────
gene_scores = ig.abs().sum(axis=1)
top20 = gene_scores.nlargest(20).index
mat   = ig.loc[top20]

# ── HIERARCHICAL CLUSTERING ─────────────────────────────
# linkage for rows (genes) and columns (latents)
row_link = linkage(mat.values, method="average", metric="correlation")
col_link = linkage(mat.values.T, method="average", metric="correlation")

# dendrogram orders
row_dendr = dendrogram(row_link, no_plot=True)
col_dendr = dendrogram(col_link, no_plot=True)
row_idx, col_idx = row_dendr["leaves"], col_dendr["leaves"]

# reorder matrix
mat_ord = mat.iloc[row_idx, :].iloc[:, col_idx]

########################################################################
# ── DRAW FIGURE ─────────────────────────────────────────
fig = plt.figure(figsize=(8, 6))
import matplotlib.gridspec as gridspec
gs = gridspec.GridSpec(2, 2,
    width_ratios=(0.15, 4),   # EVEN less space for row dendro
    height_ratios=(1, 4),
    wspace=0.0, hspace=0.0)   # NO gap

# Column dendrogram (top)
ax_col = fig.add_subplot(gs[0,1])
dendrogram(col_link, ax=ax_col, color_threshold=None,
           labels=mat_ord.columns, orientation='top')
ax_col.axis('off')

# Row dendrogram (left, trunk at right)
ax_row = fig.add_subplot(gs[1,0])
dendrogram(row_link, ax=ax_row, color_threshold=None,
           labels=mat_ord.index, orientation='left')
ax_row.axis('off')

# Heatmap
ax_heat = fig.add_subplot(gs[1,1])
im = ax_heat.matshow(mat_ord.values, aspect='auto', origin='lower',
                      cmap='viridis')

# coolwarm
# YlGnBu
# plasma
# 


# Column labels at bottom
ax_heat.set_xticks(np.arange(len(mat_ord.columns)))
ax_heat.set_xticklabels(mat_ord.columns, rotation=90, fontsize=9)
ax_heat.xaxis.set_ticks_position('bottom')
ax_heat.tick_params(axis='x', which='both',
                    labelbottom=True, labeltop=False)

# Gene labels at right
ax_heat.set_yticks(np.arange(len(mat_ord.index)))
ax_heat.set_yticklabels(mat_ord.index, fontsize=9)
ax_heat.yaxis.set_ticks_position('right')
ax_heat.tick_params(axis='y', which='both',
                    labelleft=False, labelright=True)

# Move colorbar further right to avoid overlap
cax = fig.add_axes([0.965, 0.15, 0.015, 0.7])
fig.colorbar(im, cax=cax, label='Mean IG Attribution')

plt.suptitle("Hierarchically Clustered Heatmap\nTop 20 Genes × All 50 Latents", y=0.98)
plt.tight_layout(rect=[0,0,0.96,1])

# Save
out_file = os.path.join(out_dir, "ig_top20_dendro_heatmap.pdf")
plt.savefig(out_file)
print("Saved dendrogram heatmap to", out_file)

