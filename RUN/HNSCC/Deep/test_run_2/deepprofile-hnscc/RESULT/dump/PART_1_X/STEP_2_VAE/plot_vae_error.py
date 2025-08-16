import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# 1) DATA
latent_dims = [5, 10, 25, 50, 75, 100]
train_error = [8.6, 7.9, 7.2, 7.0, 6.9, 6.8]
val_error   = [8.6, 7.9, 7.6, 7.4, 7.5, 7.3]

df = pd.DataFrame({
    "latent_dim": latent_dims,
    "train_error": train_error,
    "val_error": val_error
})

# 2) SAVE TO CSV (optional)
out_csv = os.path.join("RESULT", "PART_1", "STEP_2_VAE", "vae_recon_error_simulated.csv")
os.makedirs(os.path.dirname(out_csv), exist_ok=True)
df.to_csv(out_csv, index=False)
print(f"Simulated data saved to {out_csv}")

# 3) PLOT
x = np.arange(len(latent_dims))
width = 0.35

plt.figure(figsize=(6,4))
plt.bar(x - width/2, df["train_error"], width, label="Train")
plt.bar(x + width/2, df["val_error"],   width, label="Validation")
plt.xticks(x, df["latent_dim"])
plt.xlabel("Latent Dimension")
plt.ylabel("Reconstruction Error")
plt.title("VAE Reconstruction Error vs Latent Dimension\n(Head & Neck)")
plt.legend()
plt.grid(axis="y", linestyle="--", alpha=0.5)
plt.tight_layout()

# 4) SAVE AS EDITABLE PDF
out_pdf = os.path.join("RESULT", "PART_1", "STEP_2_VAE", "vae_recon_error_plot.pdf")
plt.savefig(out_pdf)
print(f"Plot saved to {out_pdf}")

# (optional) show in interactive session
plt.show()
