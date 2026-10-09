import pandas as pd
import os

# USER PARAMS
attr_file = "HEAD_NECK_DeepProfile_Ensemble_Gene_Importance_Weights_50L.tsv"
top_n = 20  # or 10, 50, etc.

df = pd.read_csv(attr_file, sep="\t", index_col=0)
os.makedirs("TopGenes", exist_ok=True)
csv_rows = []

for node in df.index:
    vals = df.loc[node].abs()
    top_genes = vals.sort_values(ascending=False).head(top_n).index.tolist()
    # Save as TXT (one per latent)
    with open(f"TopGenes/latent{node}_top{top_n}_genes.txt", "w") as f:
        f.write("\n".join(top_genes))
    # For CSV
    csv_rows.append([node] + top_genes)

# Save all as a CSV: each row = latent, columns = top genes
colnames = ["latent"] + [f"gene_{i+1}" for i in range(top_n)]
pd.DataFrame(csv_rows, columns=colnames).to_csv(f"TopGenes/top{top_n}_genes_per_latent.csv", index=False)
print(f"Exported top {top_n} genes for each latent as txt (TopGenes/) and csv.")
