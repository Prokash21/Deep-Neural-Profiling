import matplotlib
matplotlib.use("Agg")  # Use non-GUI backend for safe saving

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams

# === PDF + FONT SETTINGS ===
rcParams['pdf.fonttype'] = 42       # Illustrator-compatible
rcParams['font.family'] = 'Arial'   # Clean, common font

# === LOAD EXPRESSION DATA ===
df = pd.read_csv("HEAD_NECK_DATA_TOP2_JOINED_BATCH_CORRECTED_CLEANED.tsv", sep="\t", index_col=0)
n_samples, n_genes = df.shape
print(f"Loaded expression matrix with {n_samples} samples and {n_genes} genes.")

# === TRANSPOSE FOR SAMPLE-WISE BOXPLOT ===
df_t = df.T  # rows = genes, columns = samples

# === STYLISH BOXPLOT ===
plt.figure(figsize=(10, 3.5))  # Reduced width for a clean panel look

df_t.boxplot(
    showfliers=True,
    patch_artist=True,
    boxprops=dict(facecolor="#BFD3E6", color="#405396", linewidth=1),
    whiskerprops=dict(color="#4D4D4D", linewidth=1),
    capprops=dict(color="#4D4D4D", linewidth=1),
    medianprops=dict(color="#842A2A", linewidth=1.5),
    flierprops=dict(marker='o', markersize=0.01,  # tiny dots for outliers
                    markerfacecolor="#929CAF", markeredgecolor='#929CAF')
)

# Remove tick labels (sample names)
plt.xticks([], [])

# Remove grid
plt.grid(False)

# Clean labels
plt.xlabel("Samples (n=643)")
plt.ylabel("Gene (n=11020) Expression")
plt.title("Box Plot of All Samples", pad=10)

# Save high-quality vector PDF
plt.tight_layout()
plt.savefig("boxplot_all_samples_nature_style.pdf", format='pdf')
plt.close()

print("✅ Saved: boxplot_all_samples_nature_style.pdf")

