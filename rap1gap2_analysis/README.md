# RAP1GAP2 / RAP1GAP in HNSCC — multi-omic follow-up

Outputs of the cBioPortal / GDC / GEO analyses of RAP1GAP2 and its paralog RAP1GAP, used in
`zenodo_manuscript/` (Figure 8, Supplementary Tables S6–S9). Every number below comes from
the files in this folder.

## Folder layout

- `figures/` — four publication-resolution PNGs (300 dpi)
- `tables/` — small result tables suitable for main-text or supplementary tables
- `source_data/` — per-sample data behind each figure (for a "source data" supplement)
- `supplementary/` — genome-wide DESeq2 results and the raw count matrices

---

## Findings that are statistically supported

These reach significance and can be stated as positive results.

**1. Promoter methylation is inversely correlated with RAP1GAP2 mRNA (TCGA-HNSC, n = 522).**
Spearman rho = -0.26 (p = 2.07e-09); Pearson r = -0.23 on log2(RSEM+1) (p = 8.11e-08).
Direction is consistent with methylation-linked repression, but the effect is modest —
methylation explains roughly 7% of the rank variance in expression.
→ `figures/RAP1GAP2_meth_vs_expr.png`, `source_data/RAP1GAP2_methylation_expression.csv`

**2. RAP1GAP2 expression tracks CDK6 and MAPK1 positively, and MAPK3 negatively
(TCGA-HNSC, tertile split, n = 172 vs 172).**
Surviving Holm correction across the five markers tested:
CDK6 (Mann-Whitney p = 1.8e-05, Holm p = 9.2e-05, higher in RAP1GAP2-high),
MAPK1 (p = 2.3e-04, Holm p = 9.3e-04, higher in RAP1GAP2-high),
MAPK3 (p = 7.2e-04, Holm p = 2.2e-03, higher in RAP1GAP2-low).
CCND1 is nominally significant (p = 0.039) but does not survive correction
(Holm p = 0.079); CDK4 is not significant (p = 0.110).
Whole-cohort correlations agree in sign with every tertile comparison.
**Caveat for the discussion:** three of four significant markers run *opposite* to the
tumour-suppressor prediction, and MAPK1 and MAPK3 — which encode ERK2 and ERK1 — move in
opposite directions, which argues that transcript abundance is a poor proxy for pathway
activation here. Phospho-ERK staining would be needed to make a signalling claim.
→ `tables/RAP1GAP2_downstream_stats.csv`, `source_data/RAP1GAP2_downstream_correlation.csv`

**3. Copy-number alteration is common for both paralogs, but shallow (TCGA-HNSC, n = 517 with CNA calls).**
RAP1GAP2: 66 gains (12.8%) and 121 losses (23.4%).
RAP1GAP: 56 gains (10.8%) and 88 losses (17.0%).
High-level amplification is absent for both (0 and 0 samples), and deep deletion
is near-absent (1 and 0 samples). Report both denominators: broad alteration
(any gain or loss) reaches 36.3% / 29.4% of 496 fully-profiled samples, but the
strict definition (high-level amp, deep deletion or mutation) reaches only
0.4% / 1.2%. These are arm-level passenger events, not focal targeting.
→ `tables/RAP1GAP2_vs_RAP1GAP_HNSC_summary.csv`, `tables/RAP1GAP2_vs_RAP1GAP_HNSC_alterations.csv`

---

## Findings that are not significant

Each is supported by a positive control in the same analysis.

**4. Neither paralog is differentially expressed in tumour vs matched adjacent normal — in two
independent cohorts.**

| gene | TCGA-HNSC, 43 paired patients | GSE178537 (TVGH), 20 paired patients |
|---|---|---|
| RAP1GAP2 | log2FC -0.126 (95% CI -0.49 to +0.24), padj 0.57 | log2FC +0.311 (95% CI -0.08 to +0.70), padj 0.17 |
| RAP1GAP | log2FC -0.078 (95% CI -0.54 to +0.38), padj 0.79 | log2FC +0.008 (95% CI -0.71 to +0.73), padj 0.99 |

Positive controls from the same TCGA model: MMP1 and CDKN2A strongly up, KRT4 and CRNN strongly
down — all at padj < 1e-20, confirming the pipeline detects real HNSCC biology. The confidence
intervals exclude a two-fold change in either direction in both cohorts.
→ `figures/RAP1GAP_paralogs_tumor_vs_normal.png`, `figures/RAP1GAP_paralogs_LFC_two_cohorts.png`,
  `tables/RAP1GAP_paralogs_DESeq2_result.csv`, `tables/GSE178537_RAP1GAP_paralogs_DESeq2.csv`

**5. RAP1GAP2 expression does not stratify overall survival (TCGA-HNSC, tertile split).**
Log-rank p = 0.91; Cox HR (low vs high) = 1.02 (95% CI 0.74–1.40, p = 0.91).
Median OS 54.0 months (low, 78 deaths / 171) vs 49.4 months (high, 75 deaths / 172).
A continuous Cox model on the z-score was also null, so the result is not an artefact of
dichotomising. Proportional hazards was not formally assessed.
→ `figures/RAP1GAP2_survival_KM.png`, `source_data/RAP1GAP2_survival_data.csv`

**6. RAP1GAP2 is essentially never mutated in HNSCC.**
1 of 515 sequenced samples (0.19%): a single missense variant, V335A
(NM_015085.4, chr17:2898720 T>C, VAF 0.23), which falls inside the
Rap/Ran-GAP domain. RAP1GAP: 6 of 515 (1.17%). Both are far below any
plausible driver threshold.
→ `tables/RAP1GAP2_HNSC_mutations.csv`

---

## Methods

Somatic mutation, GISTIC copy-number, HM450 promoter methylation, RNA-seq expression z-scores
and clinical data for TCGA-HNSC were retrieved from the cBioPortal REST API
(studies `hnsc_tcga_pan_can_atlas_2018` and `hnsc_tcga`). Because neither cBioPortal HNSC study
contains solid-tissue normal samples, tumour-versus-normal differential expression was computed
from GDC STAR gene-level counts (TCGA-HNSC project), restricted to the 43 patients with both a
primary tumour and a matched adjacent normal. An independent replication used GSE178537
(RSEM expected counts, 20 patients with both adjacent normal and primary tumour; lymph-node
metastases excluded). Both differential-expression analyses used PyDESeq2 with the design
`~patient + condition`, retaining genes with at least 10 counts in at least 10 samples, and
Benjamini-Hochberg adjustment. Group comparisons used the Mann-Whitney U test with Holm
correction across the five marker genes; correlations are reported as Pearson (on log2 scale)
and Spearman. Survival analysis used lifelines (Kaplan-Meier, log-rank, Cox proportional
hazards).

## Data availability

- cBioPortal: https://www.cbioportal.org/api — studies `hnsc_tcga_pan_can_atlas_2018`, `hnsc_tcga`
- GDC: https://api.gdc.cancer.gov — TCGA-HNSC, "STAR - Counts"; file manifest in `source_data/gdc_paired_manifest.csv`
- GEO: GSE178537 (supplementary file `GSE178537_HNSCC_Patient_expected_count.txt.gz`)

## Limitations

- Adjacent "normal" mucosa in HNSCC may carry field-cancerisation changes, which biases a
  tumour-vs-normal contrast toward the null.
- Analyses are gene-level; isoform-level or subgroup-specific (e.g. HPV status) effects were
  not tested.
- No normal-tissue methylation is available in either cBioPortal HNSC study, so the
  tumour-vs-normal methylation comparison was not performed.
- Proportional-hazards assumptions were not formally assessed.
- GSE178537 is small for these genes (baseMean 174 and 56 vs 1342 and 1089 in TCGA); it
  replicates the TCGA null rather than independently establishing one.

---

## Update: tumour-vs-normal methylation (`methylation/` folder)

The limitation above ("no normal-tissue methylation ... tumour-versus-normal methylation
comparison was not performed") is addressed by the GDC-based analysis in
[`../methylation/`](../methylation/). It uses TCGA-HNSC HM450 data from GDC, which includes
50 solid-tissue normals alongside 528 primary tumours.

- **RAP1GAP2 promoter is hypomethylated in tumours:** 7 of 10 promoter CpGs significant,
  mean Δβ −0.20, min FDR 1.3e-47 (`methylation/results/tables/08_promoter_gene_methylation_all.csv`).
- **RAP1GAP promoter is not differentially methylated** (0 of 34 CpGs pass |Δβ| ≥ 0.2).
- Expression is unchanged (unpaired TCGA DESeq2, log2FC ≈ 0), consistent with point 4 above.

**Why the correlation differs from point 1:** point 1 reports ρ = −0.26 using tumours only,
cBioPortal's single gene-level β value and RSEM expression. The `methylation/` pipeline
reports ρ = 0.01 using tumours *and* normals, the mean of 10 promoter CpGs and log2 CPM.
Including normals (lower methylation, similar expression) flattens the within-tumour trend.
State which definition you use when reporting.
