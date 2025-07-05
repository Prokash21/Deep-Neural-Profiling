import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Load GSEA summary data
df = pd.read_csv("gsea_all_latents_summary.csv")

# Convert data types
df["FDR q-val"] = pd.to_numeric(df["FDR q-val"], errors="coerce")
df["NES"] = pd.to_numeric(df["NES"], errors="coerce")

# ── FILTER: keep only pathways with q < 0.05 in at least one latent ──
significant_terms = (
    df[df["FDR q-val"] < 0.05]["Term"]
    .dropna()
    .unique()
)
df_filtered = df[df["Term"].isin(significant_terms)]

# Pivot NES matrix
heatmap_data = df_filtered.pivot_table(index="Term", columns="latent", values="NES", aggfunc="mean")

# Annotate with p-value stars
def pval_to_stars(p):
    if pd.isna(p): return ""
    elif p < 0.0001: return "****"
    elif p < 0.001: return "***"
    elif p < 0.01: return "**"
    elif p < 0.05: return "*"
    else: return ""

annotations = df_filtered.pivot_table(index="Term", columns="latent", values="FDR q-val", aggfunc="min")
annotations = annotations.map(pval_to_stars)

# ── PLOT ──
plt.figure(figsize=(15, 0.22 * len(heatmap_data)))  # wider boxes, compact height
sns.set(style="whitegrid")

ax = sns.heatmap(
    heatmap_data,
    annot=annotations,
    fmt="",
    cmap="coolwarm",
    center=0,
    linewidths=0.4,
    linecolor='gray',
    cbar_kws={'label': 'NES'},
    square=False,
    annot_kws={"size": 8}  # slightly smaller text
)

plt.title("GSEA NES Heatmap (q < 0.05, HNSCC)", fontsize=14, pad=12)
plt.xticks(rotation=45, ha="right", fontsize=9)
plt.yticks(fontsize=7)
plt.xlabel("Latent Node", fontsize=11)
plt.ylabel("Significant Pathway", fontsize=11)
plt.tight_layout()

# Save
plt.savefig("GSEA_Significant_Only_Heatmap.pdf")
plt.savefig("GSEA_Significant_Only_Heatmap.png", dpi=300)
plt.show()
