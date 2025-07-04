import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import (roc_auc_score, average_precision_score,
                             accuracy_score, f1_score, precision_score, recall_score,
                             roc_curve, precision_recall_curve)
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout, Conv1D, MaxPooling1D, Flatten, LSTM, Input, GRU
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras.utils import to_categorical

from matplotlib import rcParams

# === Illustrator-compatible settings ===
rcParams['pdf.fonttype'] = 42   # Fonts as editable text
rcParams['font.family'] = 'Arial'

# === Load your standardized expression matrix ===
df = pd.read_csv("57_37_ready_expression_standardized.csv", index_col=0)
X = df.drop(columns=["condition"]).values
y = df["condition"].values

# === Model definitions ===
def build_mlp(input_shape):
    model = Sequential([
        Dense(64, activation='relu', input_shape=input_shape),
        Dropout(0.3),
        Dense(32, activation='relu'),
        Dense(1, activation='sigmoid')
    ])
    model.compile(optimizer=Adam(0.001), loss='binary_crossentropy', metrics=['AUC'])
    return model

def build_cnn(input_shape):
    model = Sequential([
        Input(shape=input_shape),
        Conv1D(32, 3, activation='relu'),
        MaxPooling1D(2),
        Flatten(),
        Dense(32, activation='relu'),
        Dense(1, activation='sigmoid')
    ])
    model.compile(optimizer=Adam(0.001), loss='binary_crossentropy', metrics=['AUC'])
    return model

def build_lstm(input_shape):
    model = Sequential([
        Input(shape=input_shape),
        LSTM(32),
        Dense(1, activation='sigmoid')
    ])
    model.compile(optimizer=Adam(0.001), loss='binary_crossentropy', metrics=['AUC'])
    return model

# === Wrap models for benchmarking ===
models = {
    "MLP": lambda: build_mlp((X.shape[1],)),
    "CNN": lambda: build_cnn((X.shape[1], 1)),
    "LSTM": lambda: build_lstm((X.shape[1], 1)),
}

# === Benchmarking function ===
def evaluate_model(name, model_fn, X, y, epochs=50):
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    aurocs, auprcs = [], []
    all_fpr, all_tpr, all_prec, all_rec = [], [], [], []

    for train_idx, test_idx in skf.split(X, y):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        if name in ["CNN", "LSTM"]:
            X_train = X_train[..., np.newaxis]
            X_test = X_test[..., np.newaxis]

        model = model_fn()
        model.fit(X_train, y_train, epochs=epochs, batch_size=16,
                  validation_split=0.2, verbose=0, callbacks=[EarlyStopping(patience=5)])

        y_prob = model.predict(X_test).ravel()
        y_pred = (y_prob > 0.5).astype(int)

        auroc = roc_auc_score(y_test, y_prob)
        auprc = average_precision_score(y_test, y_prob)
        aurocs.append(auroc)
        auprcs.append(auprc)

        fpr, tpr, _ = roc_curve(y_test, y_prob)
        prec, rec, _ = precision_recall_curve(y_test, y_prob)
        all_fpr.append(fpr)
        all_tpr.append(tpr)
        all_prec.append(prec)
        all_rec.append(rec)

    print(f"{name} | Mean AUROC: {np.mean(aurocs):.3f} | AUPRC: {np.mean(auprcs):.3f}")
    return {
        "model": name,
        "AUROC": np.mean(aurocs),
        "AUPRC": np.mean(auprcs),
        "fpr": all_fpr,
        "tpr": all_tpr,
        "prec": all_prec,
        "rec": all_rec
    }

# === Run all models
results = []
for name, model_fn in models.items():
    results.append(evaluate_model(name, model_fn, X, y))

# === Plotting
def plot_curves(results, kind='roc'):
    plt.figure(figsize=(8, 5.5))

    metric_key = 'AUROC' if kind == 'roc' else 'AUPRC'
    x_key = 'fpr' if kind == 'roc' else 'rec'
    y_key = 'tpr' if kind == 'roc' else 'prec'

    common_x = np.linspace(0, 1, 100)  # common base for interpolation

    for res in results:
        interpolated_curves = []

        for x_vals, y_vals in zip(res[x_key], res[y_key]):
            interp_y = np.interp(common_x, x_vals, y_vals)
            plt.plot(x_vals, y_vals, alpha=0.2)  # raw curve
            interpolated_curves.append(interp_y)

        # plot mean curve
        mean_y = np.mean(interpolated_curves, axis=0)
        label = f"{res['model']} (mean {metric_key} = {res[metric_key]:.2f})"
        plt.plot(common_x, mean_y, label=label)

    plt.title(f"Mean {metric_key} Curve")
    plt.xlabel("FPR" if kind == "roc" else "Recall")
    plt.ylabel("TPR" if kind == "roc" else "Precision")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"{kind}_comparison_plot.pdf")
    plt.show()



plot_curves(results, kind='roc')
plot_curves(results, kind='auprc')

# === Final result summary
pd.DataFrame([{k: v for k, v in r.items() if k in ['model', 'AUROC', 'AUPRC']} for r in results]).to_csv("dl_model_comparison.csv", index=False)
