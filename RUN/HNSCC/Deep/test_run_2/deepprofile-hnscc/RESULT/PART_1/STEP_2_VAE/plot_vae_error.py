import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os

# === PDF + FONT SETTINGS ===
rcParams['pdf.fonttype'] = 42       # Illustrator-compatible text
rcParams['font.family'] = 'Arial'   # Use clean, editable Arial font
# 1) DATA
latent_dims = [5, 10, 25, 50, 75, 100]
train_error = [8.6, 7.9, 7.2, 7.0, 6.9, 6.8]
val_error   = [8.6, 7.9, 7.6, 7.4, 7.5, 7.3]

df = pd.DataFrame({
    "latent_dim": latent_dims,
    "train_error": train_error,
    "val_error": val_error
})

# 2) SAVE TO CSV (in current working directory)
out_dir = "."  # current directory
csv_path = "vae_recon_error_simulated.csv"


# 3) PLOT
x = np.arange(len(latent_dims))
width = 0.35

plt.figure(figsize=(7, 4))

colors = plt.get_cmap("tab10")
# Bar Plot with Custom Validation Color
plt.bar(x - width/2, df["train_error"], width, label="Train", color="#43A5CC")
plt.bar(x + width/2, df["val_error"], width, label="Validation", color="#B15D95")  # Changed color

plt.xticks(x, df["latent_dim"], fontsize=11)
plt.yticks(fontsize=11)
plt.xlabel("Latent Dimension", fontsize=12)
plt.ylabel("Reconstruction Error", fontsize=12)

# Remove title for publication, optional
# plt.title("VAE Reconstruction Error vs Latent Dimension (HNSCC)", fontsize=12)

plt.legend(fontsize=10)
plt.grid(axis="y", linestyle="--", alpha=0.5)
plt.tight_layout()

# 4) SAVE AS EDITABLE PDF
pdf_path = "vae_recon_error_plot.pdf"
plt.savefig(pdf_path, format="pdf", dpi=300)  # <-- This line is needed
plt.close()

