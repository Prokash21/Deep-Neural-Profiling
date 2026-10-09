import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from gseapy import prerank

# ── CONFIG ─────────────────────────────────────────────────────
input_file = "HEAD_NECK_DeepProfile_Ensemble_Gene_Importance_Weights_50L.tsv"
gene_sets = "KEGG_2021_Human"  # You can also try: "GO_Biological_Process_2021"
output_dir = "GSEA_results"
os.makedirs(output_dir, exist_ok=True)

# ── LOAD ATTRIBUTION SCORES ────────────────────────────────────
df = pd.read_csv(input_file, sep="\t", index_col=0)  # rows: latent nodes, columns: genes

all_results = []
for latent in df.index:
    scores = df.loc[latent]
    scores = scores.sort_values(ascending=False)

    # Prepare ranking file for prerank
    rnk = scores.reset_index()
    rnk.columns = ["gene", "score"]
    rnk_path = os.path.join(output_dir, f"latent_{latent}_ranking.rnk")
    rnk.to_csv(rnk_path, sep="\t", index=False, header=False)

    # Run GSEA prerank
    prerank_res = prerank(
        rnk=rnk_path,
        gene_sets=gene_sets,
        outdir=None,
        min_size=5,
        max_size=1000,
        permutation_num=100,
        seed=42,
        verbose=False
    )

    res_df = prerank_res.res2d.copy()
    res_df["latent"] = latent
    all_results.append(res_df)

# ── CONCATENATE RESULTS ─────────────────────────────────────────
summary_df = pd.concat(all_results)
summary_df.reset_index(inplace=True)
summary_df.rename(columns={"Term": "name"}, inplace=True)
summary_df.to_csv(os.path.join(output_dir, "GSEA_all_latents.csv"), index=False)

# ── DOTPLOT VISUALIZATION ──────────────────────────────────────
# Ensure NES is numeric
# Convert NES to numeric
# Convert NES to numeric if it's not already
# Convert NES and FDR q-val to numeric, coerce errors to NaN
summary_df["NES"] = pd.to_numeric(summary_df["NES"], errors="coerce")
summary_df["FDR q-val"] = pd.to_numeric(summary_df["FDR q-val"], errors="coerce")

# Drop rows where FDR q-val is NaN or <= 0 to safely apply log10
summary_df = summary_df.dropna(subset=["FDR q-val"])
summary_df = summary_df[summary_df["FDR q-val"] > 0]

# Now compute -log10 of adjusted p-values
summary_df["-log10_padj"] = -np.log10(summary_df["FDR q-val"])

# Select top term per latent based on NES
top_terms = summary_df.groupby("latent").apply(lambda d: d.nlargest(1, "NES")).reset_index(drop=True)

# Plot
plt.figure(figsize=(10, max(6, len(top_terms) * 0.4)))
sns.scatterplot(
    data=top_terms,
    x="latent",
    y="name",
    size="-log10_padj",
    hue="NES",
    sizes=(50, 300),
    palette="coolwarm",
    legend="auto"
)
plt.title("Top GSEA Term per Latent (NES & FDR q-val)")
plt.xlabel("Latent Node")
plt.ylabel("GSEA Term")
plt.tight_layout()
plt.savefig("GSEA_dotplot_top_terms_per_latent.pdf")
print("Saved: GSEA_dotplot_top_terms_per_latent.pdf")
