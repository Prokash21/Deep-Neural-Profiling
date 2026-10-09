import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import rcParams
import os

# === PDF + FONT SETTINGS ===
rcParams['pdf.fonttype'] = 42       # Illustrator-compatible text
rcParams['font.family'] = 'Arial'   # Use clean, editable Arial font

# === Step 1: Load GSEA result file ===
file_path = "gsea_all_latents_summary.csv"  # Update this path if needed
df = pd.read_csv(file_path)

# === Step 2: Clean numeric columns ===
df["FDR q-val"] = pd.to_numeric(df["FDR q-val"], errors="coerce")
df["NES"] = pd.to_numeric(df["NES"], errors="coerce")

# === Step 3: Filter significant pathways (FDR q-value < 0.05) ===
significant_terms = df[df["FDR q-val"] < 0.05]["Term"].dropna().unique()
df_filtered = df[df["Term"].isin(significant_terms)]

# === Step 4: Create NES matrix ===
heatmap_data = df_filtered.pivot_table(index="Term", columns="latent", values="NES", aggfunc="mean")

# === Step 5: Select Top 50 most variable pathways ===
top50_terms = heatmap_data.var(axis=1).sort_values(ascending=False).head(50).index
top50_data = heatmap_data.loc[top50_terms]

# === Step 6: Plot ===
plt.figure(figsize=(5, 12))  # A4 landscape with a bit more height
ax = sns.heatmap(
    top50_data,
    cmap="RdBu_r",
    center=0,
    linewidths=0.5,
    linecolor='gray',
    annot=False,
    cbar_kws={"label": "NES", "shrink": 0.5, "location": "right"}
)

# === Step 7: Customize colorbar and labels ===
cbar = ax.collections[0].colorbar
cbar.ax.yaxis.set_label_position('left')
cbar.ax.yaxis.set_ticks_position('left')
cbar.set_label("NES")
cbar.ax.tick_params(labelsize=8)

# === Step 8: Titles and fonts ===
plt.title("Top 50 Enriched Pathways across Latent Variables", fontsize=12)
plt.xlabel("Latent Variable", fontsize=10)
plt.ylabel("Pathway", fontsize=10)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.tight_layout()

# === Step 9: Save as editable files ===
plt.savefig("top_50_pathways_heatmap_editable.svg")
plt.savefig("top_50_pathways_heatmap_editable.pdf")
plt.show()
