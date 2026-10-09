# utils.R
# Shared settings and helpers for 02_ and 03_ scripts.

# ---- Settings ----
# All TCGA-HNSC primary tumors (all head and neck sites) vs all solid tissue normals
FDR_CUTOFF        <- 0.05
DELTA_BETA_CUTOFF <- 0.2   # |mean tumor beta - mean normal beta|
MIN_SIG_PROBES    <- 2     # significant promoter CpGs needed to call a gene hyper/hypo
LFC_CUTOFF        <- 1     # |log2 fold change| for DEGs

PROMOTER_GROUPS <- c("TSS1500", "TSS200", "5'UTR", "1stExon")

# ---- Helpers ----
# TCGA barcode: TCGA-XX-XXXX-01A-... ; characters 14-15 are the sample type code
sample_code <- function(barcode) substr(barcode, 14, 15)
patient_id  <- function(barcode) substr(barcode, 1, 12)
sample_id   <- function(barcode) substr(barcode, 1, 15)   # e.g. TCGA-XX-XXXX-01 (ignores vial)

# Keep primary tumors (01) and solid tissue normals (11), one vial/aliquot per sample.
# Returns column indices.
select_samples <- function(barcodes) {
  which(sample_code(barcodes) %in% c("01", "11") & !duplicated(sample_id(barcodes)))
}

group_factor <- function(barcodes) {
  factor(ifelse(sample_code(barcodes) == "01", "Tumor", "Normal"),
         levels = c("Normal", "Tumor"))
}

GROUP_COLORS  <- c(Normal = "#4D9221", Tumor = "#C51B7D")
STATUS_COLORS <- c(Hyper = "#B2182B", Hypo = "#2166AC", NS = "grey70",
                   Hypermethylated = "#B2182B", Hypomethylated = "#2166AC",
                   Up = "#D6604D", Down = "#4393C3")

# ---- Packages ----
ensure_packages <- function(pkgs) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, update = FALSE)
}

# ---- Output helpers ----
# results/{figures,tables}; every figure is saved as PDF + PNG,
# every table as CSV and collected for one Excel workbook per script.
results_dir <- function() {
  d <- "results"
  for (sub in c("figures", "tables")) dir.create(file.path(d, sub), recursive = TRUE, showWarnings = FALSE)
  d
}

save_plot <- function(p, name, width = 7, height = 5) {
  for (ext in c("pdf", "png"))
    ggplot2::ggsave(file.path(results_dir(), "figures", paste0(name, ".", ext)), p,
                    width = width, height = height, dpi = 300, bg = "white")
}

# For base-graphics plots (pheatmap, UpSetR): draw_fn() is called once per device
save_base_plot <- function(draw_fn, name, width = 7, height = 5) {
  f <- file.path(results_dir(), "figures", name)
  pdf(paste0(f, ".pdf"), width = width, height = height); draw_fn(); dev.off()
  png(paste0(f, ".png"), width = width, height = height, units = "in", res = 300); draw_fn(); dev.off()
}

.excel_sheets <- list()
save_table <- function(df, name, excel = TRUE, gzip = FALSE) {
  f <- file.path(results_dir(), "tables", paste0(name, ".csv", if (gzip) ".gz"))
  con <- if (gzip) gzfile(f, "w") else f
  write.csv(df, con, row.names = FALSE)
  if (gzip) close(con)
  if (excel) .excel_sheets[[substr(name, 1, 31)]] <<- df   # Excel sheet names max 31 chars
  invisible(df)
}

write_excel <- function(file_name) {
  openxlsx::write.xlsx(.excel_sheets[order(names(.excel_sheets))], file.path(results_dir(), "tables", file_name),
                       headerStyle = openxlsx::createStyle(textDecoration = "bold"),
                       firstRow = TRUE, colWidths = "auto")
}
