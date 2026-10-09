# 01_download_tcga_hnsc.R
# Download TCGA-HNSC 450K methylation (beta values), STAR RNA-seq counts and
# clinical data from GDC, and save them as .rds files in data/.
#
# Run from the methylation/ folder:  Rscript 01_download_tcga_hnsc.R
# Disk: ~10 GB for GDCdata/ (7.6 GB methylation + ~2.4 GB RNA-seq).
# RAM: ~8 GB+ for GDCprepare on the methylation data.
# Re-running is safe: files already in GDCdata/ are skipped.

# ---- 1. Packages ----
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
pkgs <- c("TCGAbiolinks", "SummarizedExperiment", "sesame", "sesameData")   # sesame: probe annotation for GDCprepare
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p)

library(TCGAbiolinks)
library(SummarizedExperiment)

dir.create("data", showWarnings = FALSE)

# Only primary tumors and adjacent normal tissue are needed
SAMPLE_TYPES <- c("Primary Tumor", "Solid Tissue Normal")

# GDC chunk downloads sometimes arrive corrupted (tar errors) and GDCdownload can
# then stop with an error. Retry with small chunks; finished files are skipped.
download_with_retry <- function(query, files_per_chunk = 5, max_tries = 10) {
  for (i in seq_len(max_tries)) {
    ok <- tryCatch({ GDCdownload(query, files.per.chunk = files_per_chunk); TRUE },
                   error = function(e) { message("Download attempt ", i, " failed: ", conditionMessage(e)); FALSE })
    # Remove leftovers of a corrupted chunk from the working directory
    uuid_dirs <- list.files(".", pattern = "^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$")
    unlink(c(Sys.glob("*.tar.gz"), "MANIFEST.txt", uuid_dirs), recursive = TRUE)
    if (ok) return(invisible(TRUE))
    Sys.sleep(10)
  }
  stop("GDCdownload failed after ", max_tries, " attempts; re-run this script to continue.")
}

# ---- 2. DNA methylation: Illumina 450K beta values ----
query_met <- GDCquery(
  project       = "TCGA-HNSC",
  data.category = "DNA Methylation",
  data.type     = "Methylation Beta Value",
  platform      = "Illumina Human Methylation 450",
  sample.type   = SAMPLE_TYPES
)
met_samples <- getResults(query_met)
print(table(met_samples$sample_type))

download_with_retry(query_met)
met <- GDCprepare(query_met)
saveRDS(met, file.path("data", "TCGA_HNSC_methylation450_SE.rds"))

# ---- 3. RNA-seq: STAR counts (same patients, used for DEGs and correlation) ----
query_rna <- GDCquery(
  project       = "TCGA-HNSC",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type   = SAMPLE_TYPES
)
download_with_retry(query_rna)
rna <- GDCprepare(query_rna)
saveRDS(rna, file.path("data", "TCGA_HNSC_STAR_counts_SE.rds"))

# ---- 4. Clinical data (tumor site table in 02_) ----
clinical <- GDCquery_clinic(project = "TCGA-HNSC", type = "clinical")
saveRDS(clinical, file.path("data", "TCGA_HNSC_clinical.rds"))
print(table(clinical$tissue_or_organ_of_origin))

sessionInfo()
