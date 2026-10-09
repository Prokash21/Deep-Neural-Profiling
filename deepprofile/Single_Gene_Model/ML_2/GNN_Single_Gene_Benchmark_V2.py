import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve, precision_recall_curve
from sklearn.preprocessing import StandardScaler
from tensorflow.keras import layers, Model, optimizers
from spektral.layers import GCNConv

# Load and preprocess
df = pd.read_csv("GSE290057_exp_heat_20.csv", index_col=0)
features = df.drop(columns=["condition"])
labels = df["condition"].astype(int).values
scaler = StandardScaler()
features_scaled = scaler.fit_transform(features)
gene_list = features.columns

metrics_list = []
roc_data = {}
prc_data = {}

for gene in gene_list:
    X_gene = features_scaled[:, [list(gene_list).index(gene)]]  # (N, 1)
    y_all, yhat_all = [], []
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    for fold, (train_idx, test_idx) in enumerate(skf.split(X_gene, labels)):
        X_train, X_test = X_gene[train_idx], X_gene[test_idx]
        y_train, y_test = labels[train_idx], labels[test_idx]
        n_train, n_test = len(train_idx), len(test_idx)

        # Build adjacency for current split (fully connected)
        A_train = np.ones((n_train, n_train), dtype=np.float32) - np.eye(n_train, dtype=np.float32)
        A_test = np.ones((n_test, n_test), dtype=np.float32) - np.eye(n_test, dtype=np.float32)

        # Add batch dimension for full-batch mode
        X_train_batch = np.expand_dims(X_train, 0)      # (1, N_train, 1)
        A_train_batch = np.expand_dims(A_train, 0)      # (1, N_train, N_train)
        y_train_batch = np.expand_dims(y_train, 0)      # (1, N_train,)

        X_test_batch = np.expand_dims(X_test, 0)
        A_test_batch = np.expand_dims(A_test, 0)
        y_test_batch = np.expand_dims(y_test, 0)

        # Build/compile GCN for this fold (use (None, 1) and (None, None) shapes)
        N_tr = n_train
        X_in = layers.Input(shape=(N_tr, 1))
        A_in = layers.Input(shape=(N_tr, N_tr))
        out = GCNConv(8, activation="relu")([X_in, A_in])
        out = layers.Dense(1, activation="sigmoid")(out)
        model = Model(inputs=[X_in, A_in], outputs=out)
        model.compile(loss='binary_crossentropy', optimizer=optimizers.Adam(0.01))

        # Train in full batch mode (node classification)
        model.fit([X_train_batch, A_train_batch], y_train_batch, batch_size=1, epochs=50, verbose=0)

        # For test, need new model matching test shape
        N_te = n_test
        X_in_te = layers.Input(shape=(N_te, 1))
        A_in_te = layers.Input(shape=(N_te, N_te))
        out_te = GCNConv(8, activation="relu")([X_in_te, A_in_te])
        out_te = layers.Dense(1, activation="sigmoid")(out_te)
        model_te = Model(inputs=[X_in_te, A_in_te], outputs=out_te)
        model_te.compile(loss='binary_crossentropy', optimizer=optimizers.Adam(0.01))
        model_te.set_weights(model.get_weights())

        yhat = model_te.predict([X_test_batch, A_test_batch], batch_size=1)[0, :, 0]
        y_all.extend(y_test)
        yhat_all.extend(yhat)

        # Save last fold's curve
        fpr, tpr, _ = roc_curve(y_test, yhat)
        rec, prec, _ = precision_recall_curve(y_test, yhat)

    y_all = np.array(y_all)
    yhat_all = np.array(yhat_all)
    auroc = roc_auc_score(y_all, yhat_all)
    auprc = average_precision_score(y_all, yhat_all)
    metrics_list.append({'Gene': gene, 'AUROC': auroc, 'AUPRC': auprc})
    roc_data[gene] = (fpr, tpr, auroc)
    prc_data[gene] = (rec, prec, auprc)

# Save metrics
metrics_df = pd.DataFrame(metrics_list)
metrics_df_sorted_roc = metrics_df.sort_values(by="AUROC", ascending=False)
metrics_df_sorted_prc = metrics_df.sort_values(by="AUPRC", ascending=False)
metrics_df.to_csv("gnn_single_gene_metrics.csv", index=False)

# AUROC Curves
plt.figure(figsize=(10, 7))
for _, row in metrics_df_sorted_roc.iterrows():
    gene = row["Gene"]
    fpr, tpr, auc_val = roc_data[gene]
    plt.plot(fpr, tpr, label=f"{gene} (AUC={auc_val:.2f})")
plt.plot([0, 1], [0, 1], 'k--')
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("AUROC Curves - Single Gene GNN Models")
plt.legend(loc="lower right", fontsize='small')
plt.grid(True)
plt.tight_layout()
plt.savefig("gnn_single_gene_auroc_curves.pdf")
plt.close()

# AUPRC Curves
plt.figure(figsize=(10, 7))
for _, row in metrics_df_sorted_prc.iterrows():
    gene = row["Gene"]
    rec, prec, auc_pr = prc_data[gene]
    plt.plot(rec, prec, label=f"{gene} (AUPRC={auc_pr:.2f})")
plt.xlabel("Recall")
plt.ylabel("Precision")
plt.title("AUPRC Curves - Single Gene GNN Models")
plt.legend(loc="upper right", fontsize='small')
plt.grid(True)
plt.tight_layout()
plt.savefig("gnn_single_gene_auprc_curves.pdf")
plt.close()
