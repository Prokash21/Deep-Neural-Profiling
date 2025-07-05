import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib import rcParams

# === Settings for Illustrator-compatible plot ===
rcParams['pdf.fonttype'] = 42
rcParams['font.family'] = 'Arial'

# === Load g:Profiler enrichment results ===
df = pd.read_csv("all_latent_gprofiler_enrichments.csv")

# === Optional: Filter to GO:BP and KEGG (or remove this to include all sources) ===
df = df[df['source'].isin(['GO:BP', 'KEGG'])]

# === Compute -log10(p-value) for dot size and color ===
df['-log10_padj'] = -np.log10(df['p_value'])

# === Make sure 'latent' is treated as a string (categorical axis) ===
df['latent'] = df['latent'].astype(str)
df['latent'] = pd.Categorical(df['latent'], ordered=True,
                              categories=sorted(df['latent'].unique(), key=lambda x: int(x)))

# === Plot ===
plt.figure(figsize=(10, max(6, len(df['name'].unique()) * 0.25)))  # Dynamic height
sns.scatterplot(
    data=df,
    x="latent",
    y="name",
    size="-log10_padj",
    hue="-log10_padj",
    sizes=(50, 300),
    palette="viridis",
    legend="auto"
)
plt.xticks(rotation=45, ha='right')
plt.title("g:Profiler Enrichment (All Terms per Latent)")
plt.xlabel("Latent Node")
plt.ylabel("Pathway/GO Term")
plt.tight_layout()

# === Save ===
plt.savefig("all_latent_gprofiler_dotplot.pdf")
plt.savefig("all_latent_gprofiler_dotplot.png", dpi=300)
plt.savefig("all_latent_gprofiler_dotplot.svg")
print("Saved: all_latent_gprofiler_dotplot.[pdf|png|svg]")
