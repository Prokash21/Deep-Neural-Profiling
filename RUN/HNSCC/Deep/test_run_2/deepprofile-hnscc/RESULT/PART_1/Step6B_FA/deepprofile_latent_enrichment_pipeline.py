import pandas as pd
import os
from gprofiler import GProfiler
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

attr_file = "HEAD_NECK_DeepProfile_Ensemble_Gene_Importance_Weights_50L.tsv"
top_n = 20
organism = "hsapiens"

# 1. Export top genes per latent
df = pd.read_csv(attr_file, sep="\t", index_col=0)
os.makedirs("TopGenes", exist_ok=True)
gp = GProfiler(return_dataframe=True)
all_terms = []
for node in df.index:
    vals = df.loc[node].abs()
    top_genes = vals.sort_values(ascending=False).head(top_n).index.tolist()
    with open(f"TopGenes/latent{node}_top{top_n}_genes.txt", "w") as f:
        f.write("\n".join(top_genes))
    # 2. Automated g:Profiler enrichment
    enrich = gp.profile(organism=organism, query=top_genes)
    enrich['latent'] = node
    enrich = enrich.head(5)  # Keep top 5 pathways per latent
    all_terms.append(enrich)
# 3. Concatenate all enrichment results
all_terms_df = pd.concat(all_terms)
all_terms_df.to_csv("all_latent_gprofiler_enrichments.csv", index=False)
print("Saved: all_latent_gprofiler_enrichments.csv")

# 4. Prepare for dotplot: Top 5 per latent, only GO:BP and KEGG (as example)
to_plot = all_terms_df[all_terms_df['source'].isin(['GO:BP', 'KEGG'])]
to_plot['-log10_padj'] = -np.log10(to_plot['p_value'])
plt.figure(figsize=(11, max(6, len(to_plot['name'].unique()) * 0.25)))
sns.scatterplot(
    data=to_plot,
    x="latent",
    y="name",
    size="-log10_padj",
    hue="-log10_padj",
    sizes=(50, 300),
    palette="viridis",
    legend="auto"
)
plt.title(f"g:Profiler enrichment (Top {top_n} genes per latent, Top 5 terms)")
plt.xlabel("Latent Node")
plt.ylabel("Pathway/GO Term")
plt.tight_layout()
plt.savefig("latent_gprofiler_dotplot.pdf")
print("Saved: latent_gprofiler_dotplot.pdf")
