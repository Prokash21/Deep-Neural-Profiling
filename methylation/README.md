# Promoter Methylation × DEG Integration (TCGA-HNSC)

Finds genes with hyper- or hypomethylated promoters in head and neck squamous cell carcinoma (TCGA-HNSC, Illumina 450K). It then overlaps them with tumor vs normal DEGs.

## Run

Run the scripts from this folder, in order:

```bash
cd methylation
Rscript 01_download_tcga_hnsc.R        # GDC download -> data/*.rds  (~10 GB, needs ~8 GB RAM)
Rscript 02_differential_methylation.R  # DMPs + promoter-level gene calls
Rscript 03_integrate_with_degs.R       # overlap with DEGs + correlation + plots
```

Missing CRAN and Bioconductor packages are installed automatically, including ggvenn, UpSetR, pheatmap and openxlsx. All settings are in [`utils.R`](utils.R):

| Setting | Default | Meaning |
|---|---|---|
| `FDR_CUTOFF` | 0.05 | For both CpGs and DEGs |
| `DELTA_BETA_CUTOFF` | 0.2 | Minimum \|mean tumor β − mean normal β\| |
| `MIN_SIG_PROBES` | 2 | Significant promoter CpGs (all in the same direction) needed to call a gene |
| `LFC_CUTOFF` | 1 | Minimum \|log2FC\| for a DEG |

In `03_integrate_with_degs.R`, `DEG_FILE = NULL` (default) computes tumor vs normal DEGs from TCGA-HNSC RNA-seq with DESeq2, using the same samples as the methylation analysis. To use another DESeq2 table, set `DEG_FILE` to its path. TCGA RNA-seq is always used for the methylation–expression correlation, because that test needs methylation and expression from the same patients. Any DESeq2 table works if it has `log2FoldChange` and `padj` columns and gene symbols in a `symbol` or `gene` column (or the first column).

## Method

1. **Samples:** all TCGA-HNSC primary tumors (`01`, every head and neck site) and solid tissue normals (`11`), one vial per sample.
2. **Probe filtering:** CpG probes only. Probes are removed if they are on chrX/chrY, have a common SNP (MAF ≥ 0.05) at the CpG or single-base extension site, or are missing in more than 10% of samples. Any remaining missing values are filled with the probe mean.
3. **Differential methylation:** limma on M-values, tumor vs normal. Δβ is the difference in mean β.
4. **Promoter calls:** promoter probes are those annotated `TSS1500`, `TSS200`, `5'UTR` or `1stExon` (IlluminaHumanMethylation450kanno.ilmn12.hg19). A gene is called hyper- or hypomethylated if it has at least `MIN_SIG_PROBES` significant promoter CpGs in one direction and none in the other.
5. **Integration:**
   - **Hyper-Down** (possible silencing): hypermethylated and downregulated.
   - **Hypo-Up** (possible activation): hypomethylated and upregulated.
   - Each overlapping gene is also tested with a Spearman correlation between mean promoter β and log2(CPM+1) in samples that have both data types.
   - Final candidates must show the concordant pattern and a significant negative correlation.

## Outputs (`results/`)

Every figure is saved as both PDF and 300-dpi PNG in `figures/`. Every table is saved as CSV in `tables/`. Each script also writes one Excel workbook (`methylation_results.xlsx`, `integration_results.xlsx`) with one sheet per table.

**Figures**

| # | Figure | # | Figure |
|---|---|---|---|
| 01 | Beta density, per sample | 16 | Venn: hypermethylated vs downregulated |
| 02 | Beta density, by group | 17 | Venn: hypomethylated vs upregulated |
| 03 | PCA, 5,000 most variable CpGs | 18 | Venn: hypermethylated vs upregulated (discordant) |
| 04 | DMP volcano | 19 | Venn: hypomethylated vs downregulated (discordant) |
| 05 | Δβ histogram | 20 | Venn: all four sets |
| 06 | Significant CpGs by gene region | 21 | UpSet plot of the four sets |
| 07 | Significant CpGs by CpG-island context | 22 | Starburst plot (promoter Δβ vs log2FC) |
| 08 | Significant CpGs by chromosome | 23 | Category counts (Hyper-Down, Hypo-Up, ...) |
| 09 | Heatmap, top 100 DMPs | 24 | Methylation–expression correlation histogram |
| 10 | Promoter hyper/hypo gene counts | 25 | Scatter plots, top 12 candidates |
| 11 | Heatmap, top promoter genes | 26 | Heatmap, candidates: promoter β |
| 12 | Boxplots, top 20 promoter genes | 27 | Heatmap, candidates: expression z-score |
| 13 | Expression PCA (VST) | 28 | Each gene of interest: β vs expression |
| 14 | DEG volcano | 29 | Each gene of interest: tumor vs normal boxplots |
| 15 | DEG MA plot | | |

Figures 13 and 15 are only made when DEGs are computed from TCGA.

**Tables**

| # | Table |
|---|---|
| 00 | Methylation summary (sample, CpG and gene counts; cutoffs) |
| 01 | Samples by group and tissue site |
| 02 | Probe filtering summary |
| 03 | All CpGs, limma results (`.csv.gz`) |
| 04–06 | Significant, hypermethylated and hypomethylated CpGs |
| 07 | Significant CpG counts by gene region, island context and chromosome |
| 08–10 | Promoter gene-level results: all genes, hypermethylated, hypomethylated |
| 10 (integration) | Integration summary |
| 11–13 | DEGs: all genes, upregulated, downregulated |
| 14 | Venn set membership (one row per gene, TRUE/FALSE for each set) |
| 15 | Overlap enrichment (one-sided Fisher test within genes measured on both platforms) |
| 16–19 | Overlap gene lists: Hyper-Down, Hypo-Up, Hyper-Up, Hypo-Down |
| 20 | All overlapping genes, with category and correlation |
| 21–22 | Final candidates: Hyper-Down and Hypo-Up with significant negative correlation |
| 23 | Genes of interest (RAP1GAP2 by default) |

Intermediate `.rds` files (`DMP_all_probes.rds`, `promoter_beta_gene_by_sample.rds`) are kept in `results/` and used by `03_`.

## Caveats

- All head and neck sites (oral cavity, larynx, oropharynx, hypopharynx) are pooled. Table `01_sample_summary` shows how many tumors come from each site.
- HPV status strongly affects HNSCC methylation. Consider adding it as a covariate or stratifying by it.
- Cross-reactive 450K probes are not removed. If you need that, add a published list (e.g. Chen et al. 2013).
- For independent validation in oropharyngeal tissue, see GEO GSE58272: 46 OPSCC tumors and 46 adjacent normals, 27K array.
