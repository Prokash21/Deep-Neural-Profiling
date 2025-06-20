# ----------------------------------------
# Script: process_GSE178537_counts.R
# ----------------------------------------

# ---- 1. Load required libraries ----
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("R.utils", quietly = TRUE)) install.packages("R.utils")
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")

library(data.table)
library(DESeq2)
library(R.utils)

# ---- 2. Set working directory (optional) ----
# setwd("your/desired/path")  # Uncomment and set to your preferred path

# ---- 3. Download the expected count file ----
count_file_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/GSE178537/suppl/GSE178537_HNSCC_Patient_expected_count.txt.gz"
download.file(url = count_file_url, destfile = "GSE178537_expected_counts.txt.gz")

# ---- 4. Unzip the file ----
gunzip("GSE178537_expected_counts.txt.gz", overwrite = TRUE)

# ---- 5. Load the count data ----
count_data <- fread("GSE178537_expected_counts.txt")

# ---- 6. Process gene IDs ----
# Extract gene names column
gene_names <- count_data[[1]]

# Remove gene name column from data
count_data <- count_data[, -1, with = FALSE]

# Assign gene names as rownames
rownames(count_data) <- gene_names

# ---- 7. Extract Ensembl IDs and Gene Symbols ----
ensembl_ids <- sapply(strsplit(gene_names, "_"), `[`, 1)
gene_symbols <- sapply(strsplit(gene_names, "_"), `[`, 2)

# Make gene symbols unique to avoid conflicts
unique_symbols <- make.unique(gene_symbols)

# Assign unique gene symbols as row names
rownames(count_data) <- unique_symbols

# ---- 8. Combine metadata and count matrix ----
output_df <- data.frame(
  Ensembl_ID = ensembl_ids,
  Gene_Symbol = unique_symbols,
  count_data,
  check.names = FALSE
)

# ---- 9. Save processed data ----
write.csv(count_data, file = "count_data_GSE178537.csv")
write.csv(output_df, file = "count_data_with_ids.csv", row.names = FALSE)

# ---- 10. Completion message ----
message("✅ Processing complete. Files saved as:
- count_data_GSE178537.csv
- count_data_with_ids.csv")
