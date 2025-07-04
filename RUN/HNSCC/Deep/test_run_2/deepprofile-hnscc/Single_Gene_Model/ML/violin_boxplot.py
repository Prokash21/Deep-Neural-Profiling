import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams

# === Illustrator-compatible settings ===
rcParams['pdf.fonttype'] = 42   # Fonts as editable text in Illustrator
rcParams['font.family'] = 'Arial'

# === Load Data ===
df_after = pd.read_csv("57_37_ready_expression_standardized.csv", index_col=0)
df_before = pd.read_csv("57_37_merged_expression_batch.csv", index_col=0)

# Drop 'batch' column if present
df_before = df_before.drop(columns=["batch"], errors="ignore")
df_after = df_after.drop(columns=["batch"], errors="ignore")

# === Add 'Correction' label for boxplot 1 ===
df_before["Correction"] = "Before"
df_after["Correction"] = "After"

# === Combine for boxplot of batch correction effect ===
df_combined = pd.concat([df_before, df_after])
df_combined = pd.melt(df_combined, id_vars=["condition", "Correction"],
                      var_name="Gene", value_name="Expression")

# === Plot 1: Boxplot Before vs After Correction ===
plt.figure(figsize=(10, 5))
sns.boxplot(data=df_combined, x="Gene", y="Expression", hue="Correction", palette="Set1")
plt.title("Before and After Batch Effect Correction")
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig("boxplot_batch_correction_editable.pdf")
plt.close()

# === Plot 2: Boxplot by Condition (After Correction) ===
df_box = pd.melt(df_after.drop(columns=["Correction"]), id_vars="condition",
                 var_name="Gene", value_name="Expression")

plt.figure(figsize=(10, 5))
sns.boxplot(data=df_box, x="Gene", y="Expression", hue="condition", palette="Set2")
plt.title("Boxplot of Gene Expression by Condition (After Correction)")
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig("boxplot_condition_after_correction_editable.pdf")
plt.close()

# === Plot 3: Violin Plot by Condition (After Correction) ===
plt.figure(figsize=(10, 6))
sns.violinplot(
    data=df_box,
    x="Gene", y="Expression", hue="condition",
    inner="quartile",
    palette="Set2",
    scale="width",
    bw=0.2,
    cut=0
)
plt.title("Violin Plot of Gene Expression by Condition (After Correction)")
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig("violin_plot_after_batch_correction_editable.pdf")
plt.close()
