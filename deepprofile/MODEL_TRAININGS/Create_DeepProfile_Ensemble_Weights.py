###############################
# Create DeepProfile Ensemble Gene Attribution Matrices
###############################

import numpy as np
import pandas as pd
import sys
import os

# --- Inputs ---
cancer_type = sys.argv[1]        # e.g. "HEAD_NECK"
final_dim = int(sys.argv[2])     # e.g. 50 (number of ensemble latent dimensions/nodes)
num_folds = int(sys.argv[3])     # e.g. 2 (number of VAE runs you actually have)

input_folder = f'./ALL_CANCER_FILES/{cancer_type}/'
output_folder = input_folder

# --- Read VAE gene attribution matrices (per fold) ---
gene_weights_list = []
for i in range(num_folds):
    # NOTE: Use 'VAE_FILES' not 'VAE_WEIGHTS' as per your directory!
    path = os.path.join(input_folder, f'VAE_FILES/{cancer_type}_DATA_VAE_Cluster_Weights_TRAINING_{final_dim}L_fold{i}.tsv')
    print(f"Loading: {path}")
    if not os.path.exists(path):
        print(f"File does not exist: {path} — skipping.")
        continue
    data_df = pd.read_table(path, index_col=0)
    gene_weights_list.append(data_df.values.T)   # transpose to (nodes, genes)

if len(gene_weights_list) == 0:
    print("No attribution files found. Please check file paths and number of folds.")
    sys.exit(1)

# (nodes*runs, genes)
all_gene_weights = np.concatenate(gene_weights_list, axis=0)
print("Joined weights", all_gene_weights.shape)

# --- Read ensemble labels ---
labels_path = os.path.join(input_folder, f'{cancer_type}_TRAINING_DATA_kmeans_ENSEMBLE_LABELS_{final_dim}L.txt')
labels = pd.read_table(labels_path, header=None).values.flatten()
print("Ensemble labels", len(labels))

# --- Aggregate (ensemble) attribution weights by label ---
ensemble_weights = np.zeros((final_dim, all_gene_weights.shape[1]))
for label in range(final_dim):
    indices = np.where(labels == label)[0]
    if len(indices) == 0:
        print(f"Warning: No indices found for label {label}, skipping.")
        continue
    ensemble_weights[label, :] = np.mean(all_gene_weights[indices, :], axis=0)

print("Ensemble weights", ensemble_weights.shape)

# --- Save ---
gene_names = data_df.index
ensemble_df = pd.DataFrame(ensemble_weights, index=np.arange(final_dim), columns=gene_names)
out_path = os.path.join(output_folder, f'{cancer_type}_DeepProfile_Ensemble_Gene_Importance_Weights_{final_dim}L.tsv')
ensemble_df.to_csv(out_path, sep='\t')
print("Saved:", out_path)
