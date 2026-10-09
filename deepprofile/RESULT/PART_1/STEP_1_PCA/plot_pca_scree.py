import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import rcParams

# === FONT + VECTOR SETTINGS ===
rcParams['pdf.fonttype'] = 42       # Embed TrueType font for Illustrator
rcParams['font.family'] = 'Arial'   # Use Arial for all text

# === Set working directory ===
wd = os.getcwd()
print(f"Working directory: {wd}")

# === Load PCA component loadings ===
components_df = pd.read_csv(
    "HEAD_NECK_DATA_TOP2_JOINED_PCA_500L_COMPONENTS.tsv",
    sep="\t",
    index_col=0
)

# === Compute variance explained ratio ===
explained_variance = np.var(components_df.values, axis=0)
explained_variance_ratio = explained_variance / np.sum(explained_variance)

# === Full Scree Plot ===
plt.figure(figsize=(10, 5))
plt.plot(range(1, len(explained_variance_ratio) + 1),
         explained_variance_ratio, marker='o', linestyle='-')  # Default color
plt.title('Explained Variance by Principal Components (Scree Plot)')
plt.xlabel('Principal Component')
plt.ylabel('Proportion of Variance Explained')
plt.grid(True)
plt.tight_layout()
plt.savefig("pca_scree_plot.pdf", format='pdf')
plt.close()

# === Top N Scree Plot ===
N = 50
plt.figure(figsize=(8, 5))
plt.plot(range(1, N + 1),
         explained_variance_ratio[:N], marker='o', linestyle='-')  # Default color
plt.title(f'Top {N} Principal Components (Scree Plot)')
plt.xlabel('Principal Component')
plt.ylabel('Proportion of Variance Explained')
plt.grid(True)
plt.tight_layout()
plt.savefig("pca_scree_plot_top50.pdf", format='pdf')
plt.close()

print("✅ Saved: pca_scree_plot.pdf and pca_scree_plot_top50.pdf")
