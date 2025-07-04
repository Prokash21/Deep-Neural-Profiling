import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, average_precision_score
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Conv1D, Flatten, Dense, Input
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.callbacks import EarlyStopping
import matplotlib.pyplot as plt
import seaborn as sns

# === Load your data ===
df = pd.read_csv("57_37_ready_expression_standardized.csv", index_col=0)
y = df["condition"].values
gene_list = df.drop(columns=["condition"]).columns
X_all = df.drop(columns=["condition"])

# === Store results
results = []

# === Evaluate each gene
for gene in gene_list:
    X_gene = X_all[[gene]].values
    X_gene = X_gene.reshape(-1, 1, 1)  # (samples, time_steps, features)

    aurocs, auprcs = [], []
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    for train_idx, test_idx in skf.split(X_gene, y):
        X_train, X_test = X_gene[train_idx], X_gene[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        model = Sequential([
            Input(shape=(1, 1)),
            Conv1D(8, 1, activation='relu'),
            Flatten(),
            Dense(1, activation='sigmoid')
        ])
        model.compile(optimizer=Adam(0.01), loss='binary_crossentropy')

        model.fit(X_train, y_train, epochs=50, batch_size=16, verbose=0,
                  validation_split=0.2, callbacks=[EarlyStopping(patience=5, restore_best_weights=True)])

        y_prob = model.predict(X_test).ravel()
        aurocs.append(roc_auc_score(y_test, y_prob))
        auprcs.append(average_precision_score(y_test, y_prob))

    results.append({
        "Gene": gene,
        "Mean_AUROC": np.mean(aurocs),
        "Mean_AUPRC": np.mean(auprcs)
    })

# === Save and rank results
results_df = pd.DataFrame(results).sort_values(by="Mean_AUROC", ascending=False)
results_df.to_csv("cnn_single_gene_ranking.csv", index=False)

# === Plotting
plt.figure(figsize=(10, 6))
sns.barplot(data=results_df, x="Mean_AUROC", y="Gene", palette="Blues_d")
plt.title("CNN Single-Gene AUROC Ranking")
plt.tight_layout()
plt.savefig("cnn_single_gene_auroc_ranking.pdf")

plt.figure(figsize=(10, 6))
sns.barplot(data=results_df, x="Mean_AUPRC", y="Gene", palette="Greens_d")
plt.title("CNN Single-Gene AUPRC Ranking")
plt.tight_layout()
plt.savefig("cnn_single_gene_auprc_ranking.pdf")

# === Sort genes by AUROC for ranking
results_df_sorted = results_df.sort_values(by="Mean_AUROC", ascending=False).reset_index(drop=True)

# === Line plot: AUROC Ranking
plt.figure(figsize=(12, 6))
plt.plot(results_df_sorted["Gene"], results_df_sorted["Mean_AUROC"], marker='o', label="AUROC", color='blue')
plt.xticks(rotation=45, ha='right')
plt.ylabel("Mean AUROC")
plt.xlabel("Gene (sorted by AUROC)")
plt.title("CNN Single-Gene AUROC Ranking")
plt.grid(True)
plt.tight_layout()
plt.savefig("cnn_single_gene_auroc_lineplot.pdf")

# === Line plot: AUPRC Ranking
results_df_sorted_auprc = results_df.sort_values(by="Mean_AUPRC", ascending=False).reset_index(drop=True)
plt.figure(figsize=(12, 6))
plt.plot(results_df_sorted_auprc["Gene"], results_df_sorted_auprc["Mean_AUPRC"], marker='o', label="AUPRC", color='green')
plt.xticks(rotation=45, ha='right')
plt.ylabel("Mean AUPRC")
plt.xlabel("Gene (sorted by AUPRC)")
plt.title("CNN Single-Gene AUPRC Ranking")
plt.grid(True)
plt.tight_layout()
plt.savefig("cnn_single_gene_auprc_lineplot.pdf")
