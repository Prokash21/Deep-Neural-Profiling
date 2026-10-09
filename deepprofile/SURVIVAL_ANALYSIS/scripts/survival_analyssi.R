
# 1) Install / load needed package
if (!requireNamespace("survminer", quietly=TRUE)) {
  install.packages("BiocManager", repos="https://cran.rstudio.com/")
  BiocManager::install("survminer")
}
library(survival)
library(survminer)

# Read in the data
clin <- read.csv("SURVIVAL_ANALYSIS/data/clinical.csv")
sig  <- read.csv("SURVIVAL_ANALYSIS/data/signature_scores.csv")

# Merge on patient barcode
df <- merge(clin, sig, by="patient_barcode")

# Let's use the first signature as an example:
df$Sig1_high = df$Sig1_score > median(df$Sig1_score, na.rm=TRUE)

# Make a survival object
surv_obj = Surv(df$OS_time, df$OS_event)

# Plot KM curve
fit = survfit(surv_obj ~ Sig1_high, data=df)
ggsurvplot(fit, data=df, pval=TRUE,
           legend.title="Signature 1",
           legend.labs=c("Low", "High"),
           risk.table=TRUE)
