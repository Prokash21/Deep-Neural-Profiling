import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os

# === PDF + FONT SETTINGS ===
rcParams['pdf.fonttype'] = 42       # Illustrator-compatible text
rcParams['font.family'] = 'Arial'   # Use clean, editable Arial font

# === Define file paths ===
pca_file = "HEAD_NECK_DATA_TOP2_JOINED_PCA_500L.tsv"
save_file = "pca_sample_scatter_plot.pdf"

# === Load PCA-transformed data (samples × PCs) ===
pca_data_df = pd.read_csv(pca_file, sep="\t", index_col=0)
print(f"Loaded PCA data with shape: {pca_data_df.shape}")

# === Plot PC1 vs PC2 ===
plt.figure(figsize=(8, 5))
plt.scatter(
    pca_data_df.iloc[:, 0],  # PC1
    pca_data_df.iloc[:, 1],  # PC2
    alpha=0.6,
    s=18,
    edgecolor='gray',
    linewidth=0.3
)
plt.title('PCA of Samples (PC1 vs PC2)', pad=10)
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.grid(True)
plt.tight_layout()

# === Save as editable vector PDF ===
plt.savefig(save_file, format='pdf')
plt.close()

print(f"✅ Saved: {os.path.abspath(save_file)}")
