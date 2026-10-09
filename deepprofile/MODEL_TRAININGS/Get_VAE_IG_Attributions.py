import os
import numpy as np
import pandas as pd
import sys
import csv
from tensorflow.keras.models import model_from_json
from IntegratedGradients import IntegratedGradients

# Read user inputs
cancer = sys.argv[1]
dimension = int(sys.argv[2])
start = int(sys.argv[3])
end = int(sys.argv[4])

print("CANCER", cancer)
print("DIM", dimension)
print("START", start)
print("END", end)

input_folder = f'./ALL_CANCER_FILES/{cancer}/'
output_folder = f'./ALL_CANCER_FILES/{cancer}/VAE_FILES/'

# You may want to pass PCA filenames as arguments, but here's the current logic:
pca_file = input_folder + f'{cancer}_DATA_TOP2_JOINED_PCA_500L_COMPONENTS.tsv'
input_file = input_folder + f'{cancer}_DATA_TOP2_JOINED_PCA_500L.tsv'

# Load PCA weights
pca_df = pd.read_table(pca_file, index_col=0)
print("PCA COMPONENTS", pca_df.shape)
pca_components = pca_df.values

# Load input data
input_df = pd.read_table(input_file, index_col=0)
print("INPUT FILE", input_df.shape)

for vae_run in range(start, end):
    print("MODEL", vae_run)

    # Load encoder
    json_file = open(output_folder + f'VAE_{cancer}_encoder_{dimension}L_{vae_run}.json', 'r')
    loaded_model_json = json_file.read()
    json_file.close()
    encoder = model_from_json(loaded_model_json)
    encoder.load_weights(output_folder + f'VAE_{cancer}_encoder_{dimension}L_{vae_run}.weights.h5')
    print("Loaded model from disk")

    ig = IntegratedGradients(encoder)

    overall_weights = np.zeros((pca_components.shape[0], dimension))

    # For each latent node
    for latent in range(dimension):
        print("Node", latent + 1)
        weights = np.zeros((pca_components.shape[0], input_df.shape[0]))

        # For each sample
        for i, row in enumerate(input_df.values):
            ig_vals = ig.explain(row, target_neuron=latent)
            # Map PCA feature attribution back to gene space
            new_vals = np.matmul(ig_vals, pca_components.T)
            weights[:, i] = new_vals

        # Average over all samples (absolute value)
        overall_weights[:, latent] = np.mean(np.abs(weights), axis=1)

    ig_df = pd.DataFrame(overall_weights, index=pca_df.index)
    print("EXPLANATIONS DF", ig_df.shape)

    ig_df.to_csv(output_folder + f'{cancer}_DATA_VAE_Cluster_Weights_TRAINING_{dimension}L_fold{vae_run}.tsv', sep='\t', quoting=csv.QUOTE_NONE)
