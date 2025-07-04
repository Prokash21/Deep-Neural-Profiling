import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.metrics import (
    roc_auc_score, average_precision_score,
    roc_curve, precision_recall_curve
)
from scikeras.wrappers import KerasClassifier
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout
from tensorflow.keras.optimizers import Adam

# === Illustrator‐compatible settings ===
rcParams['pdf.fonttype'] = 42
rcParams['font.family'] = 'Arial'

# === Load data ===
df = pd.read_csv("57_37_ready_expression_standardized.csv", index_col=0)
X_all = df.drop(columns=["condition"])
y = df["condition"].values

# === MLP builder ===
def create_mlp_model(units=64, dropout_rate=0.3, lr=0.001):
    model = Sequential([
        Dense(units, activation='relu', input_shape=(1,)),
        Dropout(dropout_rate),
        Dense(1, activation='sigmoid')
    ])
    model.compile(optimizer=Adam(lr),
                  loss='binary_crossentropy',
                  metrics=['AUC'])
    return model

# === Prepare output folder ===
os.makedirs("SingleGenePlots", exist_ok=True)

# === Result storage ===
gene_results = []

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

for gene in X_all.columns:
    X_gene = X_all[[gene]].values

    # wrap with SciKeras
    keras_clf = KerasClassifier(
        model=create_mlp_model,
        verbose=0
    )

    # proper param grid
    param_grid = {
        'model__units':        [32, 64],
        'model__dropout_rate': [0.2, 0.3],
        'model__lr':           [0.001],
        'fit__batch_size':     [8, 16],
        'fit__epochs':         [30]
    }

    grid = GridSearchCV(
        estimator=keras_clf,
        param_grid=param_grid,
        cv=cv,
        scoring='roc_auc',
        n_jobs=-1
    )
    grid.fit(X_gene, y)

    best = grid.best_estimator_

    # collect per-fold metrics
    aurocs, auprcs = [], []
    all_fpr, all_tpr, all_prec, all_rec = [], [], [], []

    for train_idx, test_idx in cv.split(X_gene, y):
        X_tr, X_te = X_gene[train_idx], X_gene[test_idx]
        y_tr, y_te = y[train_idx], y[test_idx]

        best.fit(X_tr, y_tr)
        y_prob = best.predict_proba(X_te)[:, 1]

        aurocs.append(roc_auc_score(y_te, y_prob))
        auprcs.append(average_precision_score(y_te, y_prob))

        fpr, tpr, _ = roc_curve(y_te, y_prob)
        prec, rec, _ = precision_recall_curve(y_te, y_prob)
        all_fpr.append(fpr); all_tpr.append(tpr)
        all_prec.append(prec); all_rec.append(rec)

    # === ROC plot ===
    plt.figure(figsize=(8, 5))
    for fpr, tpr in zip(all_fpr, all_tpr):
        plt.plot(fpr, tpr, alpha=0.3)
    mean_tpr = np.mean(
        [np.interp(np.linspace(0,1,100), f, t) for f, t in zip(all_fpr, all_tpr)],
        axis=0
    )
    plt.plot(np.linspace(0,1,100), mean_tpr,
             label=f"Mean AUROC = {np.mean(aurocs):.3f}")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(f"ROC Curve – {gene}")
    plt.legend(); plt.grid(True); plt.tight_layout()
    plt.savefig(f"SingleGenePlots/{gene}_ROC.pdf")
    plt.close()

    # === PRC plot ===
    plt.figure(figsize=(8, 5))
    for rec, prec in zip(all_rec, all_prec):
        plt.plot(rec, prec, alpha=0.3)
    mean_prec = np.mean(
        [np.interp(np.linspace(0,1,100), r[::-1], p[::-1])[::-1]
         for r, p in zip(all_rec, all_prec)],
        axis=0
    )
    plt.plot(np.linspace(0,1,100), mean_prec,
             label=f"Mean AUPRC = {np.mean(auprcs):.3f}")
    plt.xlabel("Recall"); plt.ylabel("Precision")
    plt.title(f"Precision–Recall Curve – {gene}")
    plt.legend(); plt.grid(True); plt.tight_layout()
    plt.savefig(f"SingleGenePlots/{gene}_PRC.pdf")
    plt.close()

    gene_results.append({
        "Gene":        gene,
        "Mean_AUROC":  np.mean(aurocs),
        "Mean_AUPRC":  np.mean(auprcs)
    })

# === Save summary ===
pd.DataFrame(gene_results).to_csv("single_gene_mlp_results.csv", index=False)
print("Done! Single‐gene evaluation complete.")
