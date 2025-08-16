import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# Load PCA component loadings (genes × PCs)
components_df = pd.read_csv(
    "ALL_CANCER_FILES/HEAD_NECK/HEAD_NECK_DATA_TOP2_JOINED_PCA_500L_COMPONENTS.tsv",
    sep="\t",
    index_col=0
)

# Compute variance explained
explained_variance = np.var(components_df.values, axis=0)
explained_variance_ratio = explained_variance / np.sum(explained_variance)

# Scree Plot
plt.figure(figsize=(10, 5))
plt.plot(range(1, len(explained_variance_ratio) + 1), explained_variance_ratio, marker='o', linestyle='-')
plt.title('Explained Variance by Principal Components (Scree Plot)')
plt.xlabel('Principal Component')
plt.ylabel('Proportion of Variance Explained')
plt.grid(True)
plt.tight_layout()

# Save as vector PDF
save_path = r"E:\DWCT\Papia Mam\Cancer\GitHub\Cancer_Project_BMU\RUN\HNSCC\Deep\test_run_2\deepprofile-hnscc\RESULT\PART_1\STEP_1_PCA\pca_scree_plot.pdf"
plt.savefig(save_path)

# Only plot top N components
N = 50
plt.figure(figsize=(8, 6))
plt.plot(range(1, N + 1), explained_variance_ratio[:N], marker='o', linestyle='-')
plt.title(f'Top {N} Principal Components (Scree Plot)')
plt.xlabel('Principal Component')
plt.ylabel('Proportion of Variance Explained')
plt.grid(True)
plt.tight_layout()

# Save PDF in same folder
save_path = r"E:\DWCT\Papia Mam\Cancer\GitHub\Cancer_Project_BMU\RUN\HNSCC\Deep\test_run_2\deepprofile-hnscc\RESULT\PART_1\STEP_1_PCA\pca_scree_plot_top50.pdf"
plt.savefig(save_path)
