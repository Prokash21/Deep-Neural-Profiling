# ——————————————————————————
# 1) Install / load
# ——————————————————————————
if (!requireNamespace("TCGAbiolinks", quietly=TRUE))
  BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)

# ——————————————————————————
# 2) Download & prepare clinical
# ——————————————————————————
# This fetches the “patient” clinical table (all metadata)
clin <- GDCquery_clinic(project = "TCGA-HNSC", type = "patient")

# ——————————————————————————
# 3) Select & rename columns
# ——————————————————————————
# We need:
#  • patient_barcode: bcr_patient_barcode
#  • OS_time: days_to_death  (if NA, use days_to_last_follow_up)
#  • OS_event: 1 if dead, 0 if alive

# Fill in follow-up for alive patients
clin$OS_time  <- ifelse(
  is.na(clin$days_to_death),
  clin$days_to_last_follow_up,
  clin$days_to_death
)
clin$OS_event <- ifelse(clin$vital_status == "Alive", 0, 1)

# Subset & rename
out <- data.frame(
  patient_barcode = clin$bcr_patient_barcode,
  OS_time         = clin$OS_time,
  OS_event        = clin$OS_event
)

# ——————————————————————————
# 4) Write to CSV
# ——————————————————————————
dir.create("SURVIVAL_ANALYSIS/data", recursive=TRUE, showWarnings=FALSE)
write.csv(
  out,
  file = "SURVIVAL_ANALYSIS/data/clinical.csv",
  row.names = FALSE,
  quote     = FALSE
)
