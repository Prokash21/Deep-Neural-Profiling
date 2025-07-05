import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

# Load the summarized GSEA results
df = pd.read_csv("gsea_all_latents_summary.csv")

# Standardize column names
df.columns = df.columns.str.strip().str.lower()

# Strip whitespaces and convert to consistent types
df['latent'] = df['latent'].astype(str)
df['name'] = df['name'].str.strip()

# If duplicates exist, average NES scores
df_grouped = df.groupby(['name', 'latent'], as_index=False)['nes'].mean()

# Pivot to matrix: rows = pathway names, columns = latent nodes
heatmap_data = df_grouped.pivot(index='name', columns='latent', values='nes')

# Optional: sort rows by NES variability
heatmap_data = heatmap_data.loc[heatmap_data.std(axis=1).sort_values(ascending=False).index]

# Plot heatmap
plt.figure(figsize=(10, max(6, 0.35 * len(heatmap_data))))
sns.heatmap(
    heatmap_data,
    annot=True,
    fmt=".2f",
    cmap="coolwarm",
    center=0,
    cbar_kws={'label': 'Normalized Enrichment Score (NES)'},
    linewidths=0.5,
    linecolor='gray'
)

plt.title("GSEA NES Heatmap per Latent Node", fontsize=14)
plt.xlabel("Latent Node", fontsize=12)
plt.ylabel("Enriched Pathway", fontsize=12)
plt.tight_layout()

# Save output
os.makedirs("GSEA_plots", exist_ok=True)
plt.savefig("GSEA_plots/GSEA_latent_heatmap_compact.pdf")
plt.savefig("GSEA_plots/GSEA_latent_heatmap_compact.png", dpi=300)
plt.show()
