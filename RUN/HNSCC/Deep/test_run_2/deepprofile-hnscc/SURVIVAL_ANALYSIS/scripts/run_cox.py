import pandas as pd
from lifelines import CoxPHFitter
import os

clin = pd.read_csv('SURVIVAL_ANALYSIS/data/clinical.csv', index_col=0)
sig  = pd.read_csv('SURVIVAL_ANALYSIS/data/signature_scores.csv', index_col=0)
df   = clin.join(sig, how='inner')

# Prepare; rename to match lifelines’ expectations
df = df.rename(columns={'OS_time':'duration', 'OS_event':'event'})

# Collect all columns ending with _score
score_cols = [c for c in df.columns if c.lower().endswith('_score')]
cols = ['duration', 'event'] + score_cols

# Drop rows with any missing values in these columns
df_cox = df[cols].dropna()

# Fit CoxPH model
cph = CoxPHFitter()
cph.fit(df_cox, duration_col='duration', event_col='event')

# Save and print summary
os.makedirs('SURVIVAL_ANALYSIS/output/cox_tables', exist_ok=True)
summary = cph.summary
summary.to_csv('SURVIVAL_ANALYSIS/output/cox_tables/cox_summary.csv')
print(summary)


import matplotlib.pyplot as plt
import numpy as np

# Only keep signature rows (if any clinical variables were included)
summary = cph.summary
summary = summary.loc[[c for c in summary.index if c.lower().endswith('_score')]]

# Sort by HR (descending)
summary = summary.sort_values('exp(coef)', ascending=False)

# Forest plot
fig, ax = plt.subplots(figsize=(6, len(summary)*0.25 + 2))
y_pos = np.arange(len(summary))
hr = summary['exp(coef)']
err_lo = hr - summary['exp(coef) lower 95%']
err_hi = summary['exp(coef) upper 95%'] - hr

# Plot HRs with CI bars
ax.errorbar(hr, y_pos, xerr=[err_lo, err_hi], fmt='o', color='darkblue',
            ecolor='lightblue', capsize=3, lw=2)

# Reference line at HR=1
ax.axvline(1, color='grey', linestyle='--')

# Axes
ax.set_yticks(y_pos)
ax.set_yticklabels(summary.index)
ax.set_xlabel('Hazard Ratio (HR)')
ax.set_title('CoxPH: Signature Scores (multivariate)')

# Highlight significant signatures (optional)
for i, p in enumerate(summary['p']):
    if p < 0.05:
        ax.get_yticklabels()[i].set_color('red')

plt.tight_layout()
plt.gcf().subplots_adjust(left=0.30)
plt.savefig('SURVIVAL_ANALYSIS/output/cox_tables/cox_forestplot.png', dpi=300)
plt.show()
