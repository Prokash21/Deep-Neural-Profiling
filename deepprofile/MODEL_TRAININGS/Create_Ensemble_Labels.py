###############################
# This script is for learning ensemble labels for VAE models
###############################

import pandas as pd
import numpy as np
import sys
from sklearn.cluster import KMeans
import os

# ---- Read user inputs ----
cancer_type = sys.argv[1]
final_dim = int(sys.argv[2])
print("FINAL DIM " + str(final_dim))

# ---- Read all training embeddings ----
# Adjust dims and run as per your project setup
dims = [50]   # or [5, 10, 25, 50, 75, 100] if you want the full ensemble
run = 2       # or 100 for the full pipeline

data_list = []

for dim in dims:
    for i in range(run):
        print(f"Reading: dim={dim}, fold={i}")
        fname = f'../ALL_CANCER_FILES/{cancer_type}/VAE_FILES/{cancer_type}_DATA_TOP2_JOINED_encoded_{dim}L_TRAINING_fold{i}.tsv'
        data_df = pd.read_table(fname, index_col=0)
        print(data_df.shape)
        data_list.append(data_df.values)

joined_df = np.concatenate(data_list, axis=1)
print("Joined training embeddings:", joined_df.shape)

# ---- Apply kmeans clustering to samples ----
X = joined_df  # shape: (samples, features)

kmeans = KMeans(n_clusters=final_dim, random_state=123).fit(X)
print("K-means labels:", kmeans.labels_)

# ---- Save labels ----
output_dir = f'../ALL_CANCER_FILES/{cancer_type}/'
os.makedirs(output_dir, exist_ok=True)
output_file = f'{output_dir}{cancer_type}_TRAINING_DATA_kmeans_ENSEMBLE_LABELS_{final_dim}L.txt'
np.savetxt(output_file, kmeans.labels_, fmt='%d', delimiter=',')

print(f"Saved {output_file} with {len(kmeans.labels_)} labels (should match number of samples)")
