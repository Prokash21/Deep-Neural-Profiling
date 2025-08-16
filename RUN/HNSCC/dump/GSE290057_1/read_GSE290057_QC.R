


qc_data <- read.csv("GSE290057_QC.csv.gz")
head(qc_data)

write.csv(qc_data, file = "GSE290057_QC_clean.csv")
