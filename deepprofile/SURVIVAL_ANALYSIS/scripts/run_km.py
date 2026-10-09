import pandas as pd
from lifelines import KaplanMeierFitter
import matplotlib.pyplot as plt
import os

# Load data
clin = pd.read_csv('SURVIVAL_ANALYSIS/data/clinical.csv', index_col=0)
sig  = pd.read_csv('SURVIVAL_ANALYSIS/data/signature_scores.csv', index_col=0)

# Merge on patient barcode
df = clin.join(sig, how='inner')

# --- Clean NaNs ---
df = df.dropna(subset=['OS_time', 'OS_event', 'Sig1_score'])

# Ensure output directory exists
os.makedirs('SURVIVAL_ANALYSIS/output/km_plots', exist_ok=True)

# Example: split by median of signature “Sig1_score”
median = df['Sig1_score'].median()
df['group'] = df['Sig1_score'] >= median

# Fit & plot
kmf = KaplanMeierFitter()
ax = plt.subplot(111)
for label, grp in df.groupby('group'):
    kmf.fit(grp['OS_time'], event_observed=grp['OS_event'], label=f"{'High' if label else 'Low'} Sig1")
    kmf.plot_survival_function(ax=ax)
plt.title("Kaplan–Meier by Sig1_score")
plt.ylabel("Survival Probability")
plt.xlabel("Time (days)")
plt.tight_layout()
out = 'SURVIVAL_ANALYSIS/output/km_plots/sig1_km.png'
plt.savefig(out, dpi=300)
plt.close()
print("Saved KM plot to", out)
