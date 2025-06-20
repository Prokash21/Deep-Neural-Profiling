# 1. Install TCGAbiolinks if not already installed
if (!requireNamespace("TCGAbiolinks", quietly = TRUE))
  BiocManager::install("TCGAbiolinks")

library(TCGAbiolinks)

# Load the package
library(SummarizedExperiment)

# 2. Query the TCGA-HNSC RNA-seq data (STAR - Counts)
query <- GDCquery(
  project = "TCGA-HNSC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# 3. Download the data
GDCdownload(query)

# Run this in R to delete corrupted chunk
# unlink("Fri_Jun_20_23_31_14_2025_14.tar.gz")

GDCdownload(query, files.per.chunk = 20)

# 4. Prepare the expression data (creates a SummarizedExperiment object)
data <- GDCprepare(query)

# 5. Extract the counts matrix
counts <- assay(data)

# 6. Save to CSV (optional)
write.csv(counts, file = "TCGA_HNSC_STAR_Counts.csv")

# Convert metadata to a data.frame
metadata_df <- as.data.frame(colData(data))

# Identify list columns
is_list_col <- sapply(metadata_df, is.list)

# Convert list columns to character (collapse multiple values if needed)
metadata_df[is_list_col] <- lapply(metadata_df[is_list_col], function(col) {
  sapply(col, function(x) paste(as.character(x), collapse = "; "))
})

# Write to CSV
write.csv(metadata_df, "TCGA_HNSC_Metadata.csv", row.names = TRUE)
