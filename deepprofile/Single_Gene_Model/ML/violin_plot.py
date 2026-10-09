import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams

# === Illustrator-compatible settings ===
rcParams['pdf.fonttype'] = 42   # Fonts as editable text
rcParams['font.family'] = 'Arial'

# === Load batch-corrected data ===
df = pd.read_csv("57_37_ready_expression_standardized.csv", index_col=0)

# Drop batch column if it exists
df = df.drop(columns=["batch"], errors="ignore")

# Melt dataframe for plotting
df_melted = pd.melt(df, id_vars=["condition"], var_name="Gene", value_name="Expression")

# === Violin Plot ===
plt.figure(figsize=(7, 4))  # Wider for better visibility
sns.violinplot(
    data=df_melted,
    x="Gene", y="Expression",
    hue="condition",
    split=True,
    inner="quartile",
    palette="Set2",
    scale="width",
    bw=0.2
)

plt.title("Gene Expression Distribution by Condition (Standardized)")
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig("violin_plot_condition_split_editable.pdf")  # Illustrator-friendly
plt.close()
