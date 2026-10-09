# 02_differential_methylation.R
# Tumor vs normal differential methylation on TCGA-HNSC 450K data, then
# summarise to promoter-level hyper/hypomethylated genes.
# Saves all figures (PDF + PNG) and tables (CSV + one Excel workbook).
#
# Run from the methylation/ folder after 01_:  Rscript 02_differential_methylation.R
# Settings (cutoffs) are in utils.R.

source("utils.R")

# ---- 1. Packages ----
ensure_packages(c("SummarizedExperiment", "limma", "minfi",
                  "IlluminaHumanMethylation450kanno.ilmn12.hg19",
                  "ggplot2", "ggrepel", "pheatmap", "matrixStats", "openxlsx"))
library(SummarizedExperiment)
library(limma)
library(ggplot2)
library(ggrepel)

out_dir <- results_dir()
theme_set(theme_minimal(base_size = 11))

# ---- 2. Load data and select samples ----
met      <- readRDS(file.path("data", "TCGA_HNSC_methylation450_SE.rds"))
clinical <- readRDS(file.path("data", "TCGA_HNSC_clinical.rds"))

met   <- met[, select_samples(colnames(met))]
group <- group_factor(colnames(met))
print(table(group))

site <- clinical$tissue_or_organ_of_origin[match(patient_id(colnames(met)), clinical$submitter_id)]
save_table(as.data.frame(table(group = group, tissue_or_organ_of_origin = site), stringsAsFactors = FALSE) |>
             subset(Freq > 0), "01_sample_summary")

beta <- assay(met)
rm(met); invisible(gc())

# ---- 3. Probe filtering ----
suppressPackageStartupMessages(library(IlluminaHumanMethylation450kanno.ilmn12.hg19))
anno <- as.data.frame(minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))

n_start <- nrow(beta)
beta <- beta[rownames(beta) %in% rownames(anno), ]
anno <- anno[rownames(beta), ]

f_cg   <- grepl("^cg", rownames(beta))
f_sex  <- !anno$chr %in% c("chrX", "chrY")
f_snp  <- !(!is.na(anno$CpG_maf) & anno$CpG_maf >= 0.05) & !(!is.na(anno$SBE_maf) & anno$SBE_maf >= 0.05)
f_na   <- rowMeans(is.na(beta)) <= 0.10
keep   <- f_cg & f_sex & f_snp & f_na
save_table(data.frame(
  step    = c("Probes in data", "With 450K annotation", "Non-CpG (ch./rs) probes",
              "On chrX / chrY", "Common SNP at CpG / SBE (MAF >= 0.05)",
              "> 10% missing", "Probes kept"),
  n_probes = c(n_start, nrow(beta), sum(!f_cg), sum(!f_sex), sum(!f_snp), sum(!f_na), sum(keep))),
  "02_probe_filtering_summary")

beta <- beta[keep, ]
anno <- anno[keep, ]
cat("Probes kept:", nrow(beta), "of", n_start, "\n")

# Impute remaining missing values with the probe mean
na_idx <- which(is.na(beta), arr.ind = TRUE)
if (nrow(na_idx) > 0) beta[na_idx] <- rowMeans(beta, na.rm = TRUE)[na_idx[, 1]]

# ---- 4. QC figures ----
set.seed(1)
dens_probes <- sample(nrow(beta), min(20000, nrow(beta)))
dens_df <- data.frame(beta  = as.vector(beta[dens_probes, ]),
                      group = rep(group, each = length(dens_probes)),
                      sample = rep(colnames(beta), each = length(dens_probes)))
save_plot(ggplot(dens_df, aes(beta, group = sample, color = group)) +
            geom_density(linewidth = 0.2, alpha = 0.4) +
            scale_color_manual(values = GROUP_COLORS) +
            labs(x = "Beta value", y = "Density", color = NULL,
                 title = "Beta value distribution per sample"),
          "01_beta_density_per_sample")
save_plot(ggplot(dens_df, aes(beta, fill = group)) +
            geom_density(alpha = 0.5) +
            scale_fill_manual(values = GROUP_COLORS) +
            labs(x = "Beta value", y = "Density", fill = NULL, title = "Beta value distribution by group"),
          "02_beta_density_by_group")
rm(dens_df)

top_var <- order(matrixStats::rowVars(beta), decreasing = TRUE)[seq_len(min(5000, nrow(beta)))]
pca <- prcomp(t(beta[top_var, ]))
pve <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
save_plot(ggplot(data.frame(pca$x[, 1:2], group = group), aes(PC1, PC2, color = group)) +
            geom_point(size = 2, alpha = 0.8) +
            scale_color_manual(values = GROUP_COLORS) +
            labs(x = paste0("PC1 (", pve[1], "%)"), y = paste0("PC2 (", pve[2], "%)"), color = NULL,
                 title = "PCA of 5,000 most variable CpGs"),
          "03_PCA_top_variable_CpGs", 6, 5)

# ---- 5. Differential methylation (limma on M-values) ----
beta_c <- pmin(pmax(beta, 1e-3), 1 - 1e-3)
mval   <- log2(beta_c / (1 - beta_c))
rm(beta_c); invisible(gc())

design <- model.matrix(~ group)
fit    <- eBayes(lmFit(mval, design))
dmp    <- topTable(fit, coef = "groupTumor", number = Inf, sort.by = "none")
rm(mval, fit); invisible(gc())

dmp$probe       <- rownames(dmp)
dmp$mean_tumor  <- rowMeans(beta[, group == "Tumor"])
dmp$mean_normal <- rowMeans(beta[, group == "Normal"])
dmp$delta_beta  <- dmp$mean_tumor - dmp$mean_normal
dmp$status <- ifelse(dmp$adj.P.Val < FDR_CUTOFF & dmp$delta_beta >=  DELTA_BETA_CUTOFF, "Hyper",
              ifelse(dmp$adj.P.Val < FDR_CUTOFF & dmp$delta_beta <= -DELTA_BETA_CUTOFF, "Hypo", "NS"))
dmp$gene        <- anno$UCSC_RefGene_Name
dmp$gene_region <- anno$UCSC_RefGene_Group
dmp$island      <- anno$Relation_to_Island
dmp$chr         <- anno$chr
dmp$pos         <- anno$pos
dmp <- dmp[, c("probe", "gene", "gene_region", "island", "chr", "pos", "mean_normal", "mean_tumor",
               "delta_beta", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "status")]
names(dmp)[names(dmp) == "logFC"]   <- "logFC_Mvalue"
names(dmp)[names(dmp) == "AveExpr"] <- "mean_Mvalue"
print(table(dmp$status))

sig <- dmp[dmp$status != "NS", ]
sig <- sig[order(sig$adj.P.Val), ]
save_table(dmp, "03_DMP_all_probes", excel = FALSE, gzip = TRUE)
save_table(sig, "04_DMP_significant_probes")
save_table(sig[sig$status == "Hyper", ], "05_DMP_hypermethylated_probes")
save_table(sig[sig$status == "Hypo", ],  "06_DMP_hypomethylated_probes")
saveRDS(dmp, file.path(out_dir, "DMP_all_probes.rds"))

# ---- 6. DMP figures ----
set.seed(1)
plot_df <- dmp[dmp$status != "NS" | runif(nrow(dmp)) < 0.1, ]   # thin non-significant points
top_lab <- head(sig, 20)
save_plot(ggplot(plot_df, aes(delta_beta, -log10(adj.P.Val), color = status)) +
            geom_point(size = 0.6, alpha = 0.5) +
            geom_vline(xintercept = c(-DELTA_BETA_CUTOFF, DELTA_BETA_CUTOFF), linetype = "dashed") +
            geom_hline(yintercept = -log10(FDR_CUTOFF), linetype = "dashed") +
            geom_text_repel(data = top_lab, aes(label = sub(";.*", "", gene)), size = 3,
                            color = "black", max.overlaps = 30) +
            scale_color_manual(values = STATUS_COLORS,
                               labels = c(Hyper = paste0("Hyper (", sum(dmp$status == "Hyper"), ")"),
                                          Hypo  = paste0("Hypo (",  sum(dmp$status == "Hypo"),  ")"),
                                          NS    = "NS")) +
            labs(x = "Delta beta (Tumor - Normal)", y = "-log10(FDR)", color = NULL,
                 title = "Differential methylation, TCGA-HNSC tumor vs normal"),
          "04_DMP_volcano")

save_plot(ggplot(dmp, aes(delta_beta)) +
            geom_histogram(bins = 100, fill = "grey40") +
            geom_vline(xintercept = c(-DELTA_BETA_CUTOFF, DELTA_BETA_CUTOFF), linetype = "dashed", color = "red") +
            labs(x = "Delta beta (Tumor - Normal)", y = "Number of CpGs", title = "Delta beta distribution"),
          "05_delta_beta_histogram")

# Counts of significant CpGs by genomic feature
first_region <- function(x) ifelse(x == "" | is.na(x), "Intergenic", sub(";.*", "", x))
sig$region1 <- first_region(sig$gene_region)
region_tab <- as.data.frame(table(feature = sig$region1, status = sig$status))
island_tab <- as.data.frame(table(feature = sig$island,  status = sig$status))
chr_tab    <- as.data.frame(table(feature = sig$chr,     status = sig$status))
save_table(rbind(cbind(annotation = "gene_region", region_tab),
                 cbind(annotation = "CpG_island",  island_tab),
                 cbind(annotation = "chromosome",  chr_tab)), "07_DMP_counts_by_feature")

bar_plot <- function(df, title, xlab) {
  ggplot(df, aes(feature, Freq, fill = status)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = STATUS_COLORS) +
    labs(x = xlab, y = "Significant CpGs", fill = NULL, title = title) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
save_plot(bar_plot(region_tab, "Significant CpGs by gene region", NULL), "06_DMP_by_gene_region")
save_plot(bar_plot(island_tab, "Significant CpGs by CpG island context", NULL), "07_DMP_by_island_context")
chr_tab$feature <- factor(chr_tab$feature, levels = paste0("chr", 1:22))
save_plot(bar_plot(chr_tab, "Significant CpGs by chromosome", NULL), "08_DMP_by_chromosome", 9, 4.5)

# Heatmap of top 50 hyper + top 50 hypo CpGs
hm_probes <- c(head(sig$probe[sig$status == "Hyper"], 50), head(sig$probe[sig$status == "Hypo"], 50))
col_anno  <- data.frame(Group = group, row.names = colnames(beta))
ord       <- order(group)
if (length(hm_probes) >= 2) {
  save_base_plot(function() pheatmap::pheatmap(
    beta[hm_probes, ord], cluster_cols = FALSE, show_colnames = FALSE, fontsize_row = 5,
    annotation_col = col_anno, annotation_colors = list(Group = GROUP_COLORS),
    labels_row = paste(hm_probes, sub(";.*", "", dmp[hm_probes, "gene"])),
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100), breaks = seq(0, 1, length.out = 101),
    main = "Top differentially methylated CpGs (beta)", silent = FALSE),
    "09_heatmap_top_DMPs", 9, 11)
}

# ---- 7. Promoter probe -> gene map ----
# UCSC_RefGene_Name / _Group are ";"-separated and aligned (one entry per transcript)
genes   <- strsplit(anno$UCSC_RefGene_Name,  ";")
regions <- strsplit(anno$UCSC_RefGene_Group, ";")
map <- data.frame(
  probe  = rep(rownames(anno), lengths(genes)),
  gene   = unlist(genes),
  region = unlist(regions)
)
map <- unique(map[map$region %in% PROMOTER_GROUPS, c("probe", "gene")])
cat("Promoter probe-gene pairs:", nrow(map), "| genes:", length(unique(map$gene)), "\n")

# ---- 8. Gene-level promoter calls ----
pm <- merge(map, dmp[, c("probe", "delta_beta", "adj.P.Val", "status")], by = "probe")
gene_tab <- do.call(rbind, lapply(split(pm, pm$gene), function(d) data.frame(
  gene              = d$gene[1],
  n_promoter_probes = nrow(d),
  n_hyper           = sum(d$status == "Hyper"),
  n_hypo            = sum(d$status == "Hypo"),
  mean_delta_beta   = mean(d$delta_beta),
  max_abs_delta_beta= max(abs(d$delta_beta)),
  min_fdr           = min(d$adj.P.Val),
  probes_hyper      = paste(d$probe[d$status == "Hyper"], collapse = ";"),
  probes_hypo       = paste(d$probe[d$status == "Hypo"],  collapse = ";")
)))
gene_tab$promoter_status <- ifelse(gene_tab$n_hyper >= MIN_SIG_PROBES & gene_tab$n_hypo == 0, "Hypermethylated",
                            ifelse(gene_tab$n_hypo  >= MIN_SIG_PROBES & gene_tab$n_hyper == 0, "Hypomethylated", "NS"))
gene_tab <- gene_tab[order(gene_tab$promoter_status, -abs(gene_tab$mean_delta_beta)), ]
print(table(gene_tab$promoter_status))

save_table(gene_tab, "08_promoter_gene_methylation_all")
save_table(gene_tab[gene_tab$promoter_status == "Hypermethylated", ], "09_promoter_hypermethylated_genes")
save_table(gene_tab[gene_tab$promoter_status == "Hypomethylated", ],  "10_promoter_hypomethylated_genes")
save_table(data.frame(
  metric = c("Tumor samples", "Normal samples", "CpGs tested", "Hypermethylated CpGs", "Hypomethylated CpGs",
             "Genes with promoter CpGs", "Promoter-hypermethylated genes", "Promoter-hypomethylated genes",
             "FDR cutoff", "Delta beta cutoff", "Min significant promoter CpGs"),
  value  = c(sum(group == "Tumor"), sum(group == "Normal"), nrow(dmp), sum(dmp$status == "Hyper"),
             sum(dmp$status == "Hypo"), nrow(gene_tab), sum(gene_tab$promoter_status == "Hypermethylated"),
             sum(gene_tab$promoter_status == "Hypomethylated"), FDR_CUTOFF, DELTA_BETA_CUTOFF, MIN_SIG_PROBES)),
  "00_summary")

save_plot(ggplot(subset(as.data.frame(table(status = gene_tab$promoter_status)), status != "NS"),
                 aes(status, Freq, fill = status)) +
            geom_col(width = 0.6, show.legend = FALSE) +
            geom_text(aes(label = Freq), vjust = -0.3) +
            scale_fill_manual(values = STATUS_COLORS) +
            labs(x = NULL, y = "Genes", title = "Genes with differentially methylated promoters"),
          "10_promoter_gene_counts", 5, 4.5)

# ---- 9. Gene x sample mean promoter beta (used here and in 03_) ----
map <- map[map$probe %in% rownames(beta), ]
n_per_gene <- table(map$gene)
# Sum over promoter probes in chunks of samples to limit memory, then average
chunks <- split(seq_len(ncol(beta)), ceiling(seq_len(ncol(beta)) / 100))
promoter_beta <- do.call(cbind, lapply(chunks, function(cols)
  rowsum(beta[map$probe, cols, drop = FALSE], map$gene)))
promoter_beta <- promoter_beta / as.vector(n_per_gene[rownames(promoter_beta)])
saveRDS(promoter_beta, file.path(out_dir, "promoter_beta_gene_by_sample.rds"))

# Promoter-level heatmap: top 30 hyper + top 30 hypo genes
hm_genes <- c(head(gene_tab$gene[gene_tab$promoter_status == "Hypermethylated"], 30),
              head(gene_tab$gene[gene_tab$promoter_status == "Hypomethylated"], 30))
hm_genes <- intersect(hm_genes, rownames(promoter_beta))
if (length(hm_genes) >= 2) {
  save_base_plot(function() pheatmap::pheatmap(
    promoter_beta[hm_genes, ord], cluster_cols = FALSE, show_colnames = FALSE, fontsize_row = 6,
    annotation_col = col_anno, annotation_colors = list(Group = GROUP_COLORS),
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100), breaks = seq(0, 1, length.out = 101),
    main = "Top promoter hyper/hypomethylated genes (mean promoter beta)", silent = FALSE),
    "11_heatmap_top_promoter_genes", 9, 9)
}

# Boxplots: top 10 hyper + top 10 hypo genes
box_genes <- c(head(gene_tab$gene[gene_tab$promoter_status == "Hypermethylated"], 10),
               head(gene_tab$gene[gene_tab$promoter_status == "Hypomethylated"], 10))
box_genes <- intersect(box_genes, rownames(promoter_beta))
if (length(box_genes) > 0) {
  box_df <- data.frame(gene  = factor(rep(box_genes, times = ncol(promoter_beta)), levels = box_genes),
                       beta  = as.vector(promoter_beta[box_genes, , drop = FALSE]),
                       group = rep(group, each = length(box_genes)))
  save_plot(ggplot(box_df, aes(group, beta, fill = group)) +
              geom_boxplot(outlier.size = 0.4) +
              facet_wrap(~ gene, ncol = 5) +
              scale_fill_manual(values = GROUP_COLORS) +
              labs(x = NULL, y = "Mean promoter beta", fill = NULL,
                   title = "Top promoter hyper- (first rows) and hypomethylated genes"),
            "12_boxplots_top_promoter_genes", 10, 9)
}

write_excel("methylation_results.xlsx")
sessionInfo()
