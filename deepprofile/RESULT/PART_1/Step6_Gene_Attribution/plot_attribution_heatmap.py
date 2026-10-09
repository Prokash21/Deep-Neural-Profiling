import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys

# USAGE: python plot_attribution_heatmap.py <attribution_file>
# Example: python plot_attribution_heatmap.py HEAD_NECK_DeepProfile_Ensemble_Weights_50L.tsv

att_file = sys.argv[1]
df = pd.read_csv(att_file, sep="\t", index_col=0)

plt.figure(figsize=(16,8))
sns.heatmap(df, cmap="RdBu_r", center=0, cbar_kws={'label': 'Attribution Score'})
plt.xlabel("Gene")
plt.ylabel("Latent Node")
plt.title("Step 6: Latent × Gene Attribution Heatmap")
plt.tight_layout()
plt.savefig("ensemble_attribution_heatmap.pdf")
print("Saved: ensemble_attribution_heatmap.pdf")
