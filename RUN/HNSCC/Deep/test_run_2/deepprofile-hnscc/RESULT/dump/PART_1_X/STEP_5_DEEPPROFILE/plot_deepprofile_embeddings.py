import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from mpl_toolkits.mplot3d import Axes3D

latent_path = "HEAD_NECK_DeepProfile_Training_Embedding_50L.tsv"

latent = pd.read_csv(latent_path, sep="\t", index_col=0)

# Optional: load sample labels (replace with your file if you have sample categories)
# sample_labels = pd.read_csv("SAMPLE_LABELS.tsv", sep="\t", index_col=0)["label"]
# sample_labels = sample_labels.loc[latent.index]

# --- PCA 2D ---
pca = PCA(n_components=2)
latent_pca_2d = pca.fit_transform(latent)
plt.figure(figsize=(8,6))
plt.scatter(latent_pca_2d[:, 0], latent_pca_2d[:, 1], c='b', alpha=0.7)
plt.title("2D PCA of DeepProfile Ensemble Embeddings (Samples)")
plt.xlabel("PCA 1")
plt.ylabel("PCA 2")
plt.tight_layout()
plt.savefig("deepprofile_embedding_PCA_2D.png")
plt.close()

# --- PCA 3D ---
pca3d = PCA(n_components=3)
latent_pca_3d = pca3d.fit_transform(latent)
fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(latent_pca_3d[:, 0], latent_pca_3d[:, 1], latent_pca_3d[:, 2], c='b', alpha=0.7)
ax.set_title("3D PCA of DeepProfile Ensemble Embeddings (Samples)")
ax.set_xlabel("PCA 1")
ax.set_ylabel("PCA 2")
ax.set_zlabel("PCA 3")
plt.tight_layout()
plt.savefig("deepprofile_embedding_PCA_3D.png")
plt.close()

# --- t-SNE 2D ---
tsne2 = TSNE(n_components=2, random_state=42)
latent_tsne_2d = tsne2.fit_transform(latent)
plt.figure(figsize=(8,6))
plt.scatter(latent_tsne_2d[:, 0], latent_tsne_2d[:, 1], c='b', alpha=0.7)
plt.title("2D t-SNE of DeepProfile Ensemble Embeddings (Samples)")
plt.xlabel("t-SNE 1")
plt.ylabel("t-SNE 2")
plt.tight_layout()
plt.savefig("deepprofile_embedding_tSNE_2D.png")
plt.close()

# --- t-SNE 3D ---
tsne3 = TSNE(n_components=3, random_state=42)
latent_tsne_3d = tsne3.fit_transform(latent)
fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(latent_tsne_3d[:, 0], latent_tsne_3d[:, 1], latent_tsne_3d[:, 2], c='b', alpha=0.7)
ax.set_title("3D t-SNE of DeepProfile Ensemble Embeddings (Samples)")
ax.set_xlabel("t-SNE 1")
ax.set_ylabel("t-SNE 2")
ax.set_zlabel("t-SNE 3")
plt.tight_layout()
plt.savefig("deepprofile_embedding_tSNE_3D.png")
plt.close()

print("All plots saved in the current working directory.")
