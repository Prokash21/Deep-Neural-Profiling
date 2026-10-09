import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
for fold in [0, 1]:
    file = f"HEAD_NECK_DATA_VAE_Cluster_Weights_TRAINING_50L_fold{fold}.tsv"
    df = pd.read_csv(file, sep="\t", index_col=0)
    plt.figure(figsize=(11,8))
    sns.heatmap(df, cmap="vlag", center=0)
    plt.title(f"IG Attribution Heatmap Fold {fold}")
    plt.tight_layout()
    plt.savefig(f"IG_heatmap_fold{fold}.pdf")
