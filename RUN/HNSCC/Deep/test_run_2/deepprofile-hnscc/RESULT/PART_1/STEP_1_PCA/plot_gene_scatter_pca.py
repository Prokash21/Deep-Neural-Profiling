import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import numpy as np

# Try to import UMAP, else skip
try:
    import umap
    umap_available = True
except ImportError:
    umap_available = False

# === 1. LOAD DATA ===
df = pd.read_csv("HEAD_NECK_DATA_TOP2_JOINED_BATCH_CORRECTED_CLEANED.tsv", sep="\t", index_col=0)
n_samples, n_genes = df.shape
print(f"Loaded expression matrix with {n_samples} samples and {n_genes} genes.")

# === 2. Save shape as a text file ===
with open("expression_matrix_shape.txt", "w") as f:
    f.write(f"Samples: {n_samples}\nGenes: {n_genes}\n")

# === 3. Plot sample/gene count as a barplot ===
plt.figure(figsize=(8,5))
plt.bar(["Samples", "Genes"], [n_samples, n_genes], color=["steelblue", "orange"])
plt.ylabel("Count")
plt.title("Sample and Gene Count")
plt.tight_layout()
plt.savefig("matrix_sample_gene_count.pdf")
plt.close()
