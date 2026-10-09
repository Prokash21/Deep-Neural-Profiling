# plot_latent_variance.py
import pandas as pd
import matplotlib.pyplot as plt
import sys

# Load
file = sys.argv[1]
df = pd.read_csv(file, sep='\t', index_col=0)

# Variance
variances = df.var()

# Plot
plt.figure(figsize=(10, 4))
plt.bar(range(len(variances)), variances)
plt.xlabel("Latent Dimension Index")
plt.ylabel("Variance")
plt.title(f"Variance per Latent Node: {file}")
plt.tight_layout()
plt.savefig(file.replace(".tsv", "_variance.png"))
plt.savefig(file.replace(".tsv", "_variance.pdf"))
plt.show()
