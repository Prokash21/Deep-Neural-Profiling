import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.metrics import (
    roc_auc_score, average_precision_score, f1_score,
    precision_score, recall_score, accuracy_score,
    roc_curve, precision_recall_curve
)
from sklearn.preprocessing import StandardScaler
import umap.umap_ as umap
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Conv1D, Flatten, Dense, Input
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.callbacks import EarlyStopping
from scikeras.wrappers import KerasClassifier
from matplotlib import rcParams

# === Illustrator-friendly settings ===
rcParams['pdf.fonttype'] = 42
rcParams['font.family'] = 'Arial'

# === Load Data ===
df = pd.read_csv("57_37_ready_expression_standardized.csv", index_col=0)
y = df["condition"].values
gene_list = df.drop(columns=["condition"]).columns
X_all = df.drop(columns=["condition"])

# === Preprocessing: Scale gene expression ===
features = df.drop(columns=["condition"])
scaled_features = StandardScaler().fit_transform(features)
df_scaled = pd.DataFrame(scaled_features, columns=features.columns, index=df.index)
df_scaled["condition"] = df["condition"].values

# === UMAP visualization ===
reducer = umap.UMAP(random_state=42)
embedding = reducer.fit_transform(scaled_features)
plt.figure(figsize=(6, 5))
sns.scatterplot(x=embedding[:, 0], y=embedding[:, 1], hue=df["condition"], palette="coolwarm", s=50)
plt.title("UMAP of Standardized Expression (Genes)")
plt.xlabel("UMAP 1")
plt.ylabel("UMAP 2")
plt.legend(title="Condition", labels=["Control (0)", "HNSCC (1)"])
plt.tight_layout()
plt.savefig("umap_20_gene_expression.pdf")
plt.close()

# === Use scaled data for modeling ===
y = df_scaled["condition"].values
gene_list = df_scaled.drop(columns=["condition"]).columns
X_all = df_scaled.drop(columns=["condition"])

# === Output storage ===
gene_metrics = []
roc_data = {}
prc_data = {}

# === CNN model builder ===
def create_cnn_model(filters=8, learning_rate=0.01, **kwargs):
    model = Sequential([
        Input(shape=(1, 1)),
        Conv1D(filters=filters, kernel_size=1, activation='relu'),
        Flatten(),
        Dense(1, activation='sigmoid')
    ])
    model.compile(
        optimizer=Adam(learning_rate=learning_rate),
        loss="binary_crossentropy",
        metrics=["accuracy"]
    )
    return model

# === Loop through each gene ===
for gene in gene_list:
    X_gene = X_all[[gene]].values.reshape(-1, 1, 1)

    model = KerasClassifier(
        model=create_cnn_model,
        verbose=0
    )

    param_grid = {
        "model__filters": [4, 8, 16],
        "model__learning_rate": [0.001, 0.01],
        "batch_size": [16],
        "epochs": [50]
    }

    grid = GridSearchCV(estimator=model, param_grid=param_grid, cv=3, scoring='roc_auc', n_jobs=-1)
    grid_result = grid.fit(X_gene, y)
    best_model = grid_result.best_estimator_

    # Cross-validation for evaluation
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    y_true_all, y_prob_all = [], []

    for train_idx, test_idx in skf.split(X_gene, y):
        X_train, X_test = X_gene[train_idx], X_gene[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        best_model.fit(X_train, y_train, epochs=50, batch_size=16, verbose=0,
                       callbacks=[EarlyStopping(patience=5, restore_best_weights=True)],
                       validation_split=0.2)

        y_prob = best_model.predict(X_test)
        if y_prob.ndim > 1:
            y_prob = y_prob[:, 0]  # For binary case, get the probability of class 1

        y_true_all.extend(y_test)
        y_prob_all.extend(y_prob)

    y_pred = (np.array(y_prob_all) > 0.5).astype(int)

    auroc = roc_auc_score(y_true_all, y_prob_all)
    auprc = average_precision_score(y_true_all, y_prob_all)
    acc = accuracy_score(y_true_all, y_pred)
    f1 = f1_score(y_true_all, y_pred)
    prec = precision_score(y_true_all, y_pred)
    rec = recall_score(y_true_all, y_pred)

    gene_metrics.append({
        "Gene": gene,
        "AUROC": auroc,
        "AUPRC": auprc,
        "Accuracy": acc,
        "F1": f1,
        "Precision": prec,
        "Recall": rec
    })

    fpr, tpr, _ = roc_curve(y_true_all, y_prob_all)
    precision, recall, _ = precision_recall_curve(y_true_all, y_prob_all)
    roc_data[gene] = (fpr, tpr, auroc)
    prc_data[gene] = (recall, precision, auprc)

# === Metrics Summary ===
metrics_df = pd.DataFrame(gene_metrics)
metrics_df_sorted_roc = metrics_df.sort_values(by="AUROC", ascending=False)
metrics_df_sorted_prc = metrics_df.sort_values(by="AUPRC", ascending=False)
metrics_df.to_csv("cnn_single_gene_metrics.csv", index=False)

# === AUROC Curve Plot ===
plt.figure(figsize=(8, 5))
for _, row in metrics_df_sorted_roc.iterrows():
    gene = row["Gene"]
    fpr, tpr, auc_val = roc_data[gene]
    plt.plot(fpr, tpr, label=f"{gene} (AUC={auc_val:.2f})")
plt.plot([0, 1], [0, 1], 'k--')
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("AUROC Curves - Single Gene CNN Models")
plt.legend(loc="lower right", fontsize='small')
plt.grid(True)
plt.tight_layout()
plt.savefig("cnn_single_gene_auroc_curves.pdf")
plt.close()

# === AUPRC Curve Plot ===
plt.figure(figsize=(8, 5))
for _, row in metrics_df_sorted_prc.iterrows():
    gene = row["Gene"]
    rec, prec, auc_pr = prc_data[gene]
    plt.plot(rec, prec, label=f"{gene} (AUPRC={auc_pr:.2f})")
plt.xlabel("Recall")
plt.ylabel("Precision")
plt.title("AUPRC Curves - Single Gene CNN Models")
plt.legend(loc="upper right", fontsize='small')
plt.grid(True)
plt.tight_layout()
plt.savefig("cnn_single_gene_auprc_curves.pdf")
plt.close()
