# 1) Install / load needed package
if (!requireNamespace("TCGAbiolinks", quietly=TRUE)) {
  install.packages("BiocManager", repos="https://cran.rstudio.com/")
  BiocManager::install("TCGAbiolinks")
}
library(TCGAbiolinks)

# 2) Download & prepare clinical
clin <- GDCquery_clinic(project = "TCGA-HNSC", type = "clinical")

# 3) Compute OS time & event
clin$OS_time  <- ifelse(
  is.na(clin$days_to_death),
  clin$days_to_last_follow_up,
  clin$days_to_death
)
clin$OS_event <- ifelse(clin$vital_status == "Alive", 0, 1)

out <- data.frame(
  patient_barcode = clin$bcr_patient_barcode,
  OS_time         = clin$OS_time,
  OS_event        = clin$OS_event
)

# 4) Write CSV
dir.create("SURVIVAL_ANALYSIS/data", recursive=TRUE, showWarnings=FALSE)
write.csv(
  out,
  file = "SURVIVAL_ANALYSIS/data/clinical.csv",
  row.names = FALSE,
  quote     = FALSE
)
message("▶ clinical.csv written to SURVIVAL_ANALYSIS/data/")
