import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import matplotlib.pyplot as plt
import os

# --- INPUT FILES ---
CLINICAL = 'SURVIVAL_ANALYSIS/data/clinical.csv'
SCORES   = 'SURVIVAL_ANALYSIS/data/signature_scores.csv'
OUTPUT   = 'SURVIVAL_ANALYSIS/output/km_plots_all/'
os.makedirs(OUTPUT, exist_ok=True)

# --- LOAD DATA ---
clin = pd.read_csv(CLINICAL, index_col=0)
sig  = pd.read_csv(SCORES, index_col=0)
df = clin.join(sig, how='inner')

# --- FIND SIGNATURE SCORE COLUMNS ---
score_cols = [c for c in df.columns if c.lower().endswith('_score') or c.isdigit()]

# --- COLLECT SUMMARY RESULTS ---
summary = []

for sig_col in score_cols:
    # Median split
    median = df[sig_col].median()
    df['group'] = df[sig_col] >= median

    # Ensure types are numeric
    df['OS_time'] = pd.to_numeric(df['OS_time'], errors='coerce')
    df['OS_event'] = pd.to_numeric(df['OS_event'], errors='coerce')

    # Filter groups for plotting and stats
    kmf = KaplanMeierFitter()
    ax = plt.figure(figsize=(7,6)).gca()
    for label, grp in df.groupby('group'):
        name = f"{'High' if label else 'Low'} {sig_col}"
        grp = grp.dropna(subset=['OS_time', 'OS_event'])  # <---- KEY FIX
        if grp.empty:
            print(f"Group {name} is empty after dropping NaNs, skipping...")
            continue
        kmf.fit(grp['OS_time'], event_observed=grp['OS_event'], label=name)
        kmf.plot_survival_function(ax=ax)
    
    # Prepare for log-rank test (also drop NaNs)
    high = df[df['group']].dropna(subset=['OS_time', 'OS_event'])
    low  = df[~df['group']].dropna(subset=['OS_time', 'OS_event'])
    if len(high) == 0 or len(low) == 0:
        pval = float('nan')
        print(f"Skipping log-rank for {sig_col}, one group is empty after dropping NaNs.")
    else:
        res  = logrank_test(high['OS_time'], low['OS_time'],
                            event_observed_A=high['OS_event'],
                            event_observed_B=low['OS_event'])
        pval = res.p_value

    # Save plot (PNG + text-editable PDF)
    plt.title(f"KM: {sig_col} (p={pval:.3g})")
    plt.ylabel("Survival Probability")
    plt.xlabel("Time (days)")
    plt.legend()
    plt.tight_layout()
    png_path = os.path.join(OUTPUT, f'{sig_col}_km.png')
    pdf_path = os.path.join(OUTPUT, f'{sig_col}_km.pdf')
    plt.savefig(png_path, dpi=300)
    plt.savefig(pdf_path, format='pdf')  # PDF = vector, editable text
    plt.close()
    
    print(f"Saved {sig_col}: p={pval:.3g}")
    summary.append({'signature': sig_col, 'p_value': pval,
                    'high_median': high['OS_time'].median(),
                    'low_median': low['OS_time'].median(),
                    'n_high': len(high), 'n_low': len(low)})


# --- SUMMARY TABLE ---
summary_df = pd.DataFrame(summary)
summary_csv = os.path.join(OUTPUT, 'km_signature_summary.csv')
summary_df.to_csv(summary_csv, index=False)
print(f"\nSaved summary to {summary_csv}")
