import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.metrics import (
    accuracy_score, roc_auc_score, average_precision_score, f1_score,
    precision_score, recall_score, confusion_matrix, roc_curve,
    precision_recall_curve
)
from sklearn.neural_network import MLPClassifier
from sklearn.inspection import permutation_importance
from sklearn.base import clone
from matplotlib import rcParams

# === Illustrator-compatible settings ===
rcParams['pdf.fonttype'] = 42   # Fonts as editable text
rcParams['font.family'] = 'Arial'

# STEP 1: Load preprocessed data
df = pd.read_csv("57_37_ready_expression_standardized.csv", index_col=0)
X = df.drop(columns=['condition'])
y = df['condition']

# STEP 2: Set up model and hyperparameter grid
import seaborn as sns  # <- Add this
from sklearn.neural_network import MLPClassifier

# Optionally adjust the model:
mlp = MLPClassifier(
    hidden_layer_sizes=(100,),
    max_iter=2000,        # increase if needed
    early_stopping=True,  # stop automatically if no improvement
    random_state=42
)
param_grid = {
    'hidden_layer_sizes': [(50,), (100,), (100, 50)],
    'alpha': [0.0001, 0.001],
    'learning_rate_init': [0.001, 0.01]
}
cv_outer = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

metrics_summary = []
mean_fpr = np.linspace(0, 1, 100)
tpr_list, prc_list = [], []

# STEP 3: Cross-validation loop with hyperparameter tuning
for train_idx, test_idx in cv_outer.split(X, y):
    X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
    y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]

    clf = GridSearchCV(mlp, param_grid, cv=3, scoring='roc_auc', n_jobs=-1)
    clf.fit(X_train, y_train)
    best_model = clf.best_estimator_

    y_pred = best_model.predict(X_test)
    y_prob = best_model.predict_proba(X_test)[:, 1]

    # Scores
    metrics_summary.append({
        "Accuracy": accuracy_score(y_test, y_pred),
        "AUROC": roc_auc_score(y_test, y_prob),
        "AUPRC": average_precision_score(y_test, y_prob),
        "F1": f1_score(y_test, y_pred),
        "Precision": precision_score(y_test, y_pred),
        "Recall": recall_score(y_test, y_pred)
    })

    # ROC and PRC curves
    fpr, tpr, _ = roc_curve(y_test, y_prob)
    tpr_list.append(np.interp(mean_fpr, fpr, tpr))

    precision, recall, _ = precision_recall_curve(y_test, y_prob)
    pr_interp = np.interp(mean_fpr, recall[::-1], precision[::-1])
    prc_list.append(pr_interp)

# STEP 4: Plot AUROC
plt.figure()
for tpr in tpr_list:
    plt.plot(mean_fpr, tpr, color='lightblue', alpha=0.3)
plt.plot(mean_fpr, np.mean(tpr_list, axis=0), label='Mean ROC', color='blue')
plt.plot([0, 1], [0, 1], 'k--')
plt.xlabel("FPR")
plt.ylabel("TPR")
plt.title("AUROC Curve")
plt.legend()
plt.savefig("auroc_curve_nn.pdf")

# Plot AUPRC
plt.figure()
for pr in prc_list:
    plt.plot(mean_fpr, pr, color='lightgreen', alpha=0.3)
plt.plot(mean_fpr, np.mean(prc_list, axis=0), label='Mean PRC', color='green')
plt.xlabel("Recall")
plt.ylabel("Precision")
plt.title("AUPRC Curve")
plt.legend()
plt.savefig("auprc_curve_nn.pdf")

# STEP 5: Finalize model on all data
final_model = GridSearchCV(mlp, param_grid, cv=5, scoring='roc_auc', n_jobs=-1)
final_model.fit(X, y)
best_final_model = final_model.best_estimator_

# STEP 6: Feature Importance (Permutation)
result = permutation_importance(best_final_model, X, y, scoring='roc_auc', n_repeats=30, random_state=42)
importances = pd.DataFrame({
    'Feature': X.columns,
    'Importance': result.importances_mean
}).sort_values(by='Importance', ascending=False)

# Plot importance
plt.figure(figsize=(8, 5))
sns.barplot(data=importances, y='Feature', x='Importance')
plt.title("Permutation Feature Importance (AUROC)")
plt.tight_layout()
plt.savefig("nn_feature_importance.pdf")

# STEP 7: Single-Gene Evaluation
single_gene_results = []
for gene in X.columns:
    X_single = X[[gene]]
    clf_single = clone(best_final_model)
    scores = []

    for train_idx, test_idx in StratifiedKFold(n_splits=5, shuffle=True, random_state=42).split(X_single, y):
        X_tr, X_te = X_single.iloc[train_idx], X_single.iloc[test_idx]
        y_tr, y_te = y.iloc[train_idx], y.iloc[test_idx]
        clf_single.fit(X_tr, y_tr)
        y_prob_single = clf_single.predict_proba(X_te)[:, 1]
        scores.append(roc_auc_score(y_te, y_prob_single))
    
    single_gene_results.append({
        'Gene': gene,
        'Mean AUROC': np.mean(scores),
        'Std AUROC': np.std(scores)
    })

# Save outputs
pd.DataFrame(metrics_summary).to_csv("nn_cv_metrics_summary.csv", index=False)
importances.to_csv("nn_feature_importance.csv", index=False)
pd.DataFrame(single_gene_results).to_csv("nn_single_gene_evaluation.csv", index=False)
