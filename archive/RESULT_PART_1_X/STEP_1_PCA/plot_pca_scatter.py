import pandas as pd
import matplotlib.pyplot as plt

# Use absolute path to PCA-transformed data (samples × PCs)
pca_file = r"E:\DWCT\Papia Mam\Cancer\GitHub\Cancer_Project_BMU\RUN\HNSCC\Deep\test_run_2\deepprofile-hnscc\ALL_CANCER_FILES\HEAD_NECK\HEAD_NECK_DATA_TOP2_JOINED_PCA_500L.tsv"

# Load PCA data
pca_data_df = pd.read_csv(pca_file, sep="\t", index_col=0)

# Scatter plot of PC1 vs PC2
plt.figure(figsize=(8, 6))
plt.scatter(pca_data_df.iloc[:, 0], pca_data_df.iloc[:, 1], alpha=0.6)
plt.title('PCA of Samples (PC1 vs PC2)')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.grid(True)
plt.tight_layout()

# Save plot as PDF
save_path = r"E:\DWCT\Papia Mam\Cancer\GitHub\Cancer_Project_BMU\RUN\HNSCC\Deep\test_run_2\deepprofile-hnscc\RESULT\PART_1\STEP_1_PCA\pca_sample_scatter_plot.pdf"
plt.savefig(save_path)
