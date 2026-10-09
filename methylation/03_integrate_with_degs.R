# 03_integrate_with_degs.R
# Overlap promoter hyper/hypomethylated genes (from 02_) with DEGs and check that
# promoter methylation is negatively correlated with expression in matched samples.
# Saves all figures (PDF + PNG), Venn/UpSet diagrams and tables (CSV + one Excel workbook).
#
# Run from the methylation/ folder after 02_:  Rscript 03_integrate_with_degs.R

source("utils.R")

# ---- 1. Packages ----
ensure_packages(c("SummarizedExperiment", "DESeq2", "ggplot2", "ggrepel", "pheatmap",
                  "ggvenn", "UpSetR", "openxlsx"))
library(SummarizedExperiment)
library(DESeq2)
library(ggplot2)
library(ggrepel)

out_dir <- results_dir()
theme_set(theme_minimal(base_size = 11))

# DEG source:
#   NULL      -> compute TCGA-HNSC tumor vs normal DEGs with DESeq2 (same samples as 02_)
#   "file.csv"-> use your own DESeq2 table with log2FoldChange and padj columns, and gene
#                symbols in a "symbol" or "gene" column (otherwise the first column is used).
# TCGA RNA-seq is always used for the methylation-expression correlation (same patients).
DEG_FILE <- NULL

GENES_OF_INTEREST <- c("RAP1GAP2")

# ---- 2. Load methylation results and TCGA expression ----
gene_tab      <- read.csv(file.path(out_dir, "tables", "08_promoter_gene_methylation_all.csv"))
promoter_beta <- readRDS(file.path(out_dir, "promoter_beta_gene_by_sample.rds"))

rna      <- readRDS(file.path("data", "TCGA_HNSC_STAR_counts_SE.rds"))
rna      <- rna[, select_samples(colnames(rna))]
group    <- group_factor(colnames(rna))
print(table(group))

counts <- assay(rna, "unstranded")
symbol <- rowData(rna)$gene_name

log_cpm <- log2(t(t(counts) / colSums(counts)) * 1e6 + 1)
rownames(log_cpm) <- symbol
log_cpm <- log_cpm[!duplicated(symbol) & !is.na(symbol), ]

# ---- 3. DEGs ----
if (is.null(DEG_FILE)) {
  keep <- rowSums(counts >= 10) >= min(table(group))
  dds  <- DESeqDataSetFromMatrix(counts[keep, ], data.frame(group = group), design = ~ group)
  dds  <- DESeq(dds)
  res  <- as.data.frame(results(dds, contrast = c("group", "Tumor", "Normal")))
  res$gene_id <- rownames(res)
  res$gene    <- symbol[keep]

  # Expression PCA (VST)
  vsd <- assay(if (nrow(dds) >= 1000) vst(dds, blind = TRUE) else varianceStabilizingTransformation(dds, blind = TRUE))
  top <- order(matrixStats::rowVars(vsd), decreasing = TRUE)[seq_len(min(2000, nrow(vsd)))]
  pca <- prcomp(t(vsd[top, ]))
  pve <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  save_plot(ggplot(data.frame(pca$x[, 1:2], group = group), aes(PC1, PC2, color = group)) +
              geom_point(size = 2, alpha = 0.8) +
              scale_color_manual(values = GROUP_COLORS) +
              labs(x = paste0("PC1 (", pve[1], "%)"), y = paste0("PC2 (", pve[2], "%)"), color = NULL,
                   title = "Expression PCA (VST, 2,000 most variable genes)"),
            "13_expression_PCA", 6, 5)
  rm(vsd)
} else {
  res <- read.csv(DEG_FILE, check.names = FALSE)
  sym_col  <- intersect(c("symbol", "gene", "gene_name", "Gene"), names(res))[1]
  res$gene <- if (!is.na(sym_col)) res[[sym_col]] else res[[1]]
}

# One row per gene symbol: keep the most significant entry
res <- res[!is.na(res$padj) & !is.na(res$gene) & res$gene != "", ]
res <- res[order(res$padj), ]
res <- res[!duplicated(res$gene), ]
res$deg_status <- ifelse(res$padj < FDR_CUTOFF & res$log2FoldChange >=  LFC_CUTOFF, "Up",
                  ifelse(res$padj < FDR_CUTOFF & res$log2FoldChange <= -LFC_CUTOFF, "Down", "NS"))
print(table(res$deg_status))

res <- res[, c("gene", setdiff(names(res), "gene"))]
save_table(res, "11_DEG_all_genes", excel = FALSE)
save_table(res[res$deg_status == "Up", ],   "12_DEG_upregulated")
save_table(res[res$deg_status == "Down", ], "13_DEG_downregulated")

# DEG volcano and MA plots
deg_lab <- rbind(head(res[res$deg_status != "NS", ], 20), res[res$gene %in% GENES_OF_INTEREST, ])
deg_lab <- deg_lab[!duplicated(deg_lab$gene), ]
deg_colors <- c(Up = "#D6604D", Down = "#4393C3", NS = "grey70")
save_plot(ggplot(res, aes(log2FoldChange, -log10(padj), color = deg_status)) +
            geom_point(size = 0.6, alpha = 0.5) +
            geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = "dashed") +
            geom_hline(yintercept = -log10(FDR_CUTOFF), linetype = "dashed") +
            geom_text_repel(data = deg_lab, aes(label = gene), size = 3, color = "black", max.overlaps = 30) +
            scale_color_manual(values = deg_colors,
                               labels = c(Up   = paste0("Up (",   sum(res$deg_status == "Up"),   ")"),
                                          Down = paste0("Down (", sum(res$deg_status == "Down"), ")"),
                                          NS   = "NS")) +
            labs(x = "log2 fold change (Tumor / Normal)", y = "-log10(FDR)", color = NULL,
                 title = "Differential expression"),
          "14_DEG_volcano")
if ("baseMean" %in% names(res)) {
  save_plot(ggplot(res, aes(log10(baseMean + 1), log2FoldChange, color = deg_status)) +
              geom_point(size = 0.5, alpha = 0.5) +
              geom_hline(yintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = "dashed") +
              scale_color_manual(values = deg_colors) +
              labs(x = "log10(mean normalized count + 1)", y = "log2 fold change", color = NULL,
                   title = "MA plot"),
            "15_DEG_MA_plot")
}

# ---- 4. Overlap ----
merged <- merge(gene_tab, res[, c("gene", "log2FoldChange", "padj", "deg_status")], by = "gene")
merged$category <- with(merged, ifelse(
  promoter_status == "Hypermethylated" & deg_status == "Down", "Hyper-Down",
  ifelse(promoter_status == "Hypomethylated" & deg_status == "Up", "Hypo-Up",
  ifelse(promoter_status == "Hypermethylated" & deg_status == "Up", "Hyper-Up",
  ifelse(promoter_status == "Hypomethylated" & deg_status == "Down", "Hypo-Down", "Other")))))
print(table(merged$category))

sets <- list(
  Hypermethylated = gene_tab$gene[gene_tab$promoter_status == "Hypermethylated"],
  Hypomethylated  = gene_tab$gene[gene_tab$promoter_status == "Hypomethylated"],
  Upregulated     = res$gene[res$deg_status == "Up"],
  Downregulated   = res$gene[res$deg_status == "Down"])

# Gene x set membership table (one row per gene in any set)
all_genes  <- sort(unique(unlist(sets)))
membership <- data.frame(gene = all_genes, sapply(sets, function(s) all_genes %in% s))
save_table(membership, "14_venn_set_membership")

# Fisher test for each overlap, within genes measured on both platforms
universe <- intersect(gene_tab$gene, res$gene)
fisher_row <- function(a, b) {
  A <- intersect(sets[[a]], universe); B <- intersect(sets[[b]], universe)
  k <- length(intersect(A, B))
  tab <- matrix(c(k, length(A) - k, length(B) - k, length(universe) - length(union(A, B))), 2)
  ft <- fisher.test(tab, alternative = "greater")
  data.frame(set_1 = a, set_2 = b, universe = length(universe), n_set_1 = length(A), n_set_2 = length(B),
             overlap = k, expected = round(length(A) * length(B) / length(universe), 2),
             odds_ratio = unname(ft$estimate), p_value = ft$p.value)
}
enrich <- rbind(fisher_row("Hypermethylated", "Downregulated"), fisher_row("Hypomethylated", "Upregulated"),
                fisher_row("Hypermethylated", "Upregulated"),   fisher_row("Hypomethylated", "Downregulated"))
enrich$fdr <- p.adjust(enrich$p_value, method = "BH")
save_table(enrich, "15_overlap_enrichment_fisher")

# ---- 5. Venn and UpSet diagrams ----
venn <- function(s, title, fill) {
  ggvenn::ggvenn(s, fill_color = fill, stroke_size = 0.4, set_name_size = 4, text_size = 4) +
    labs(title = title) + theme(plot.title = element_text(hjust = 0.5))
}
save_plot(venn(sets[c("Hypermethylated", "Downregulated")],
               "Promoter hypermethylated vs downregulated", c("#B2182B", "#4393C3")),
          "16_venn_hyper_vs_down", 6, 5)
save_plot(venn(sets[c("Hypomethylated", "Upregulated")],
               "Promoter hypomethylated vs upregulated", c("#2166AC", "#D6604D")),
          "17_venn_hypo_vs_up", 6, 5)
save_plot(venn(sets[c("Hypermethylated", "Upregulated")],
               "Promoter hypermethylated vs upregulated (discordant)", c("#B2182B", "#D6604D")),
          "18_venn_hyper_vs_up", 6, 5)
save_plot(venn(sets[c("Hypomethylated", "Downregulated")],
               "Promoter hypomethylated vs downregulated (discordant)", c("#2166AC", "#4393C3")),
          "19_venn_hypo_vs_down", 6, 5)
save_plot(venn(sets, "Promoter methylation and DEG sets",
               c("#B2182B", "#2166AC", "#D6604D", "#4393C3")),
          "20_venn_all_four_sets", 8, 7)

if (sum(lengths(sets) > 0) >= 2) {
  save_base_plot(function() print(UpSetR::upset(
    UpSetR::fromList(sets[lengths(sets) > 0]), order.by = "freq", nsets = 4,
    sets.bar.color = "grey45", main.bar.color = "grey30", text.scale = 1.3)),
    "21_upset_methylation_DEG_sets", 8, 5)
}

# ---- 6. Methylation-expression correlation in matched samples ----
met_ids <- sample_id(colnames(promoter_beta))
rna_ids <- sample_id(colnames(log_cpm))
shared  <- intersect(met_ids[!duplicated(met_ids)], rna_ids[!duplicated(rna_ids)])
cat("Samples with both methylation and RNA-seq:", length(shared), "\n")
shared_group <- group_factor(colnames(promoter_beta)[match(shared, met_ids)])

test_genes <- intersect(
  union(merged$gene[merged$category != "Other"], GENES_OF_INTEREST),
  intersect(rownames(promoter_beta), rownames(log_cpm)))

b <- promoter_beta[test_genes, match(shared, met_ids), drop = FALSE]
e <- log_cpm[test_genes, match(shared, rna_ids), drop = FALSE]
cor_tab <- do.call(rbind, lapply(test_genes, function(g) {
  ct <- suppressWarnings(cor.test(b[g, ], e[g, ], method = "spearman"))
  data.frame(gene = g, spearman_rho = unname(ct$estimate), cor_p = ct$p.value)
}))
if (is.null(cor_tab)) cor_tab <- data.frame(gene = character(), spearman_rho = numeric(), cor_p = numeric())
cor_tab$cor_fdr <- p.adjust(cor_tab$cor_p, method = "BH")

merged <- merge(merged, cor_tab, by = "gene", all.x = TRUE)
merged <- merged[order(merged$category, merged$cor_fdr), ]
merged <- merged[, c("gene", "category", "promoter_status", "deg_status", "mean_delta_beta",
                     "log2FoldChange", "padj", "spearman_rho", "cor_fdr",
                     setdiff(names(merged), c("gene", "category", "promoter_status", "deg_status",
                                              "mean_delta_beta", "log2FoldChange", "padj",
                                              "spearman_rho", "cor_fdr")))]

for (cat_name in c("Hyper-Down", "Hypo-Up", "Hyper-Up", "Hypo-Down"))
  save_table(merged[merged$category == cat_name, ],
             paste0(c("Hyper-Down" = "16", "Hypo-Up" = "17", "Hyper-Up" = "18", "Hypo-Down" = "19")[cat_name],
                    "_overlap_", gsub("-", "_", cat_name)))
save_table(merged, "20_overlap_all_genes", excel = FALSE)

# Strongest candidates: concordant direction and significant negative correlation
candidates <- merged[merged$category %in% c("Hyper-Down", "Hypo-Up") &
                     !is.na(merged$cor_fdr) & merged$cor_fdr < FDR_CUTOFF & merged$spearman_rho < 0, ]
save_table(candidates[candidates$category == "Hyper-Down", ], "21_candidates_hyper_down")
save_table(candidates[candidates$category == "Hypo-Up", ],    "22_candidates_hypo_up")
cat("Candidates (concordant + negative correlation):\n"); print(table(candidates$category))

save_table(data.frame(
  metric = c("Samples (methylation & RNA)", "Genes on both platforms", "Promoter-hypermethylated genes",
             "Promoter-hypomethylated genes", "Upregulated DEGs", "Downregulated DEGs",
             "Hyper-Down", "Hypo-Up", "Hyper-Up (discordant)", "Hypo-Down (discordant)",
             "Hyper-Down candidates (rho < 0, FDR < cutoff)", "Hypo-Up candidates (rho < 0, FDR < cutoff)",
             "DEG source"),
  value = c(length(shared), length(universe), lengths(sets),
            sum(merged$category == "Hyper-Down"), sum(merged$category == "Hypo-Up"),
            sum(merged$category == "Hyper-Up"),   sum(merged$category == "Hypo-Down"),
            sum(candidates$category == "Hyper-Down"), sum(candidates$category == "Hypo-Up"),
            if (is.null(DEG_FILE)) "TCGA DESeq2" else DEG_FILE)),
  "10_integration_summary")

# ---- 7. Integration figures ----
cat_colors <- c("Hyper-Down" = "#B2182B", "Hypo-Up" = "#2166AC",
                "Hyper-Up" = "#F4A582", "Hypo-Down" = "#92C5DE", Other = "grey75")

lab <- rbind(head(candidates[order(candidates$cor_fdr), ], 20), merged[merged$gene %in% GENES_OF_INTEREST, ])
lab <- lab[!duplicated(lab$gene), ]
save_plot(ggplot(merged, aes(mean_delta_beta, log2FoldChange, color = category)) +
            geom_point(size = 0.8, alpha = 0.6) +
            geom_vline(xintercept = c(-DELTA_BETA_CUTOFF, DELTA_BETA_CUTOFF), linetype = "dashed") +
            geom_hline(yintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = "dashed") +
            geom_text_repel(data = lab, aes(label = gene), size = 3, color = "black", max.overlaps = 30) +
            scale_color_manual(values = cat_colors) +
            labs(x = "Mean promoter delta beta (Tumor - Normal)", y = "log2 fold change (Tumor / Normal)",
                 color = NULL, title = "Starburst: promoter methylation vs expression (TCGA-HNSC)"),
          "22_starburst_methylation_vs_expression", 7, 6)

cat_counts <- as.data.frame(table(category = factor(merged$category, levels = names(cat_colors))))
save_plot(ggplot(subset(cat_counts, category != "Other"), aes(category, Freq, fill = category)) +
            geom_col(width = 0.6, show.legend = FALSE) +
            geom_text(aes(label = Freq), vjust = -0.3) +
            scale_fill_manual(values = cat_colors) +
            labs(x = NULL, y = "Genes", title = "Promoter methylation / expression categories"),
          "23_category_counts", 6, 4.5)

if (nrow(cor_tab) > 0) {
  cor_plot_df <- merge(cor_tab, merged[, c("gene", "category")], by = "gene")
  save_plot(ggplot(cor_plot_df, aes(spearman_rho, fill = category)) +
              geom_histogram(bins = 40, position = "identity", alpha = 0.7) +
              geom_vline(xintercept = 0, linetype = "dashed") +
              scale_fill_manual(values = cat_colors) +
              labs(x = "Spearman rho (promoter beta vs expression)", y = "Genes", fill = NULL,
                   title = "Methylation-expression correlation of overlapping genes"),
            "24_correlation_histogram")
}

# Scatter plots of the top 12 candidates
top_c <- head(candidates$gene[order(candidates$cor_fdr)], 12)
if (length(top_c) > 0) {
  sc_df <- data.frame(gene  = factor(rep(top_c, times = length(shared)), levels = top_c),
                      beta  = as.vector(b[top_c, , drop = FALSE]),
                      expr  = as.vector(e[top_c, , drop = FALSE]),
                      group = rep(shared_group, each = length(top_c)))
  save_plot(ggplot(sc_df, aes(beta, expr)) +
              geom_point(aes(color = group), size = 0.8, alpha = 0.6) +
              geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5, formula = y ~ x) +
              facet_wrap(~ gene, scales = "free", ncol = 4) +
              scale_color_manual(values = GROUP_COLORS) +
              labs(x = "Mean promoter beta", y = "log2(CPM + 1)", color = NULL,
                   title = "Top candidates: promoter methylation vs expression"),
            "25_scatter_top_candidates", 11, 8)
}

# Paired heatmaps of candidates: promoter methylation and expression (z-score), same sample order
hm_genes <- head(candidates$gene[order(candidates$category, candidates$cor_fdr)], 60)
if (length(hm_genes) >= 2) {
  ord      <- order(shared_group)
  col_anno <- data.frame(Group = shared_group, row.names = shared)
  row_anno <- data.frame(Category = candidates$category[match(hm_genes, candidates$gene)], row.names = hm_genes)
  ann_cols <- list(Group = GROUP_COLORS, Category = cat_colors[c("Hyper-Down", "Hypo-Up")])
  bm <- b[hm_genes, ord, drop = FALSE]; colnames(bm) <- shared[ord]
  em <- e[hm_genes, ord, drop = FALSE]; colnames(em) <- shared[ord]
  em <- t(scale(t(em)))
  save_base_plot(function() pheatmap::pheatmap(
    bm, cluster_cols = FALSE, cluster_rows = FALSE, show_colnames = FALSE, fontsize_row = 6,
    annotation_col = col_anno, annotation_row = row_anno, annotation_colors = ann_cols,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100), breaks = seq(0, 1, length.out = 101),
    main = "Candidates: mean promoter beta", silent = FALSE),
    "26_heatmap_candidates_methylation", 9, 9)
  save_base_plot(function() pheatmap::pheatmap(
    pmin(pmax(em, -3), 3), cluster_cols = FALSE, cluster_rows = FALSE, show_colnames = FALSE, fontsize_row = 6,
    annotation_col = col_anno, annotation_row = row_anno, annotation_colors = ann_cols,
    color = colorRampPalette(c("#4393C3", "white", "#D6604D"))(100), breaks = seq(-3, 3, length.out = 101),
    main = "Candidates: expression (z-score of log2 CPM)", silent = FALSE),
    "27_heatmap_candidates_expression", 9, 9)
}

# ---- 8. Genes of interest ----
goi_rows <- list()
for (g in GENES_OF_INTEREST) {
  cat("\n====", g, "====\n")
  goi_rows[[g]] <- data.frame(
    gene = g,
    promoter_status = c(gene_tab$promoter_status[gene_tab$gene == g], NA)[1],
    mean_delta_beta = c(gene_tab$mean_delta_beta[gene_tab$gene == g], NA)[1],
    n_promoter_probes = c(gene_tab$n_promoter_probes[gene_tab$gene == g], NA)[1],
    log2FoldChange  = c(res$log2FoldChange[res$gene == g], NA)[1],
    padj            = c(res$padj[res$gene == g], NA)[1],
    deg_status      = c(res$deg_status[res$gene == g], NA)[1],
    spearman_rho    = c(cor_tab$spearman_rho[cor_tab$gene == g], NA)[1],
    cor_fdr         = c(cor_tab$cor_fdr[cor_tab$gene == g], NA)[1])
  print(goi_rows[[g]])

  if (g %in% rownames(b)) {
    df <- data.frame(beta = b[g, ], expr = e[g, ], group = shared_group)
    save_plot(ggplot(df, aes(beta, expr, color = group)) +
                geom_point(alpha = 0.6) +
                geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black", formula = y ~ x) +
                scale_color_manual(values = GROUP_COLORS) +
                labs(x = paste(g, "mean promoter beta"), y = paste(g, "log2(CPM + 1)"), color = NULL,
                     title = paste0(g, ": rho = ", round(goi_rows[[g]]$spearman_rho, 2))),
              paste0("28_", g, "_methylation_vs_expression"), 5.5, 4.5)
    long <- rbind(data.frame(value = df$beta, group = df$group, measure = "Mean promoter beta"),
                  data.frame(value = df$expr, group = df$group, measure = "log2(CPM + 1)"))
    save_plot(ggplot(long, aes(group, value, fill = group)) +
                geom_boxplot(outlier.size = 0.5, show.legend = FALSE) +
                facet_wrap(~ measure, scales = "free_y") +
                scale_fill_manual(values = GROUP_COLORS) +
                labs(x = NULL, y = NULL, title = paste(g, "in tumor vs normal")),
              paste0("29_", g, "_boxplots"), 7, 4)
  }
}
save_table(do.call(rbind, goi_rows), "23_genes_of_interest")

write_excel("integration_results.xlsx")
sessionInfo()
