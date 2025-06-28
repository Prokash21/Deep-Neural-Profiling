import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# ─── Configuration ──────────────────────────────────────────────
base = r"E:\DWCT\Papia Mam\Cancer\GitHub\Cancer_Project_BMU\RUN\HNSCC\Deep\test_run_2\deepprofile-hnscc"
ig_folder = os.path.join(base, "ALL_CANCER_FILES", "HEAD_NECK", "VAE_FILES")
out_dir   = os.path.join(base, "RESULT", "PART_1", "STEP_3_IG")
os.makedirs(out_dir, exist_ok=True)

# ─── Load & Average Folds ──────────────────────────────────────
paths = [
    os.path.join(ig_folder, "HEAD_NECK_DATA_VAE_Cluster_Weights_TRAINING_50L_fold0.tsv"),
    os.path.join(ig_folder, "HEAD_NECK_DATA_VAE_Cluster_Weights_TRAINING_50L_fold1.tsv"),
]
dfs = [pd.read_csv(p, sep="\t", index_col=0) for p in paths]
# ensure same shape & columns
assert dfs[0].shape == dfs[1].shape
assert all(dfs[0].columns == dfs[1].columns)

avg_df = (dfs[0] + dfs[1]) / 2.0

# ─── 1) Heatmap: Top 20 genes × first 10 features ────────────
# pick top 20 genes by overall attribution magnitude
gene_scores = avg_df.abs().sum(axis=1)
top_genes   = gene_scores.nlargest(20).index
dims        = avg_df.columns[:10]

heat_data = avg_df.loc[top_genes, dims].values

plt.figure(figsize=(8, 6))
im = plt.imshow(heat_data, aspect="auto", interpolation="nearest")
plt.colorbar(im, label="Mean IG Attribution")
plt.xticks(np.arange(len(dims)), dims, rotation=45)
plt.yticks(np.arange(len(top_genes)), top_genes)
plt.title("Top 20 Genes × First 10 Latent Features\n(Avg. IG Attribution)")
plt.tight_layout()
plt.savefig(os.path.join(out_dir, "ig_heatmap_top20x10.pdf"))

# ─── 2) Bar-plot: Top 15 genes for the first latent feature ───
first_dim = avg_df.columns[0]
bar_data  = avg_df[first_dim].abs().nlargest(15)

plt.figure(figsize=(6, 5))
y_pos = np.arange(len(bar_data))
plt.barh(y_pos, bar_data.values)
plt.yticks(y_pos, bar_data.index)
plt.gca().invert_yaxis()
plt.xlabel("Mean |IG Attribution|")
plt.title(f"Top 15 Genes for {first_dim}")
plt.tight_layout()
plt.savefig(os.path.join(out_dir, f"ig_bar_top15_{first_dim}.pdf"))

print("Step 3 plots saved to:", out_dir)
