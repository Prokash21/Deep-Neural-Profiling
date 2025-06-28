import pandas as pd
import matplotlib.pyplot as plt
import umap
import os

base = r"E:\DWCT\Papia Mam\Cancer\GitHub\Cancer_Project_BMU\RUN\HNSCC\Deep\test_run_2\deepprofile-hnscc"

# ← Corrected paths to the encoded latent files
enc0 = os.path.join(
    base,
    "ALL_CANCER_FILES",
    "HEAD_NECK",
    "VAE_FILES",
    "HEAD_NECK_DATA_TOP2_JOINED_encoded_50L_TRAINING_fold0.tsv",
)
enc1 = os.path.join(
    base,
    "ALL_CANCER_FILES",
    "HEAD_NECK",
    "VAE_FILES",
    "HEAD_NECK_DATA_TOP2_JOINED_encoded_50L_TRAINING_fold1.tsv",
)

# Load and label
df0 = pd.read_csv(enc0, sep="\t", index_col=0)
df0["fold"] = "fold0"
df1 = pd.read_csv(enc1, sep="\t", index_col=0)
df1["fold"] = "fold1"
combined = pd.concat([df0, df1])

# UMAP reduction
reducer = umap.UMAP(n_components=2, random_state=42)
embed = reducer.fit_transform(combined.iloc[:, :-1])

# Plot
plt.figure(figsize=(8, 6))
for fold, color in zip(["fold0", "fold1"], ["C0", "C1"]):
    mask = combined["fold"] == fold
    plt.scatter(embed[mask, 0], embed[mask, 1], label=fold, alpha=0.6, s=10)
plt.legend(title="Fold")
plt.title("UMAP of VAE 50-D Latent Space")
plt.xlabel("UMAP1")
plt.ylabel("UMAP2")
plt.grid(True)
plt.tight_layout()

# Save
out = os.path.join(base, "RESULT", "PART_1", "STEP_2_VAE", "vae_latent_umap.pdf")
plt.savefig(out)
