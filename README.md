# Deep Neural Profiling Reveals RAP1GAP2 as a Latent Regulator in Head and Neck Cancer

Code, analyses and manuscript for a project that began as my **4th-year undergraduate project (June–July 2025)** in the
Department of Biochemistry and Molecular Biology, Shahjalal University of Science and Technology, and was extended with
a TCGA-HNSC multi-omic follow-up. The deep learning workflow follows the **DeepProfile** framework
(Qiu et al., *Nature Biomedical Engineering* 2025, <https://doi.org/10.1038/s41551-024-01290-8>).

**Authors:** Joy Prokash Debnath, Tanvir Hossain, Papia Rahman

**Preprint:** the manuscript is in [`zenodo_manuscript/`](zenodo_manuscript/) (Zenodo DOI to be added after upload).

---

## Abstract

**Background:** Conventional differential gene expression (DGE) analysis inadequately captures the complex molecular
changes that drive the progression of head and neck cancer, including oropharyngeal carcinoma. Variational
autoencoders (VAEs) combined with Integrated Gradients, as in the DeepProfile framework, yield interpretable latent
spaces from large cancer expression compendia.

**Methods:** Following DeepProfile, the head and neck cancer expression compendium (643 samples from 26 GEO datasets,
11,020 genes) was compressed to 500 principal components and used to train a VAE with 50 latent variables. Integrated
Gradients was used to determine the contribution of each gene to each latent variable. Genes with consistently high
attribution were selected as candidate regulators and characterised by pathway enrichment, GSEA, supervised deep
learning classifiers and DGE analysis in independent RNA-seq data. The top candidate, RAP1GAP2, was further examined
across DNA methylation, copy number, mutation, paired tumour–normal expression and survival in TCGA-HNSC and
GSE178537, together with a genome-wide integration of promoter methylation and gene expression.

**Results:** RAP1GAP2 was among the 20 genes with the highest mean attribution and was the most important feature in
the supervised classifiers: an MLP trained on the 18 measured candidate genes reached a mean AUPRC of 0.86 and AUROC of
0.80, and a RAP1GAP2-only classifier reached a mean AUPRC of 0.769. This occurred despite the lack of differential
expression in tumours relative to matched normal tissue (TCGA-HNSC, 43 pairs; GSE178537, 20 pairs). The RAP1GAP2
promoter was hypomethylated in tumours (mean Δβ −0.20), but methylation was only weakly associated with expression
(Spearman ρ = −0.26). RAP1GAP2 was rarely mutated, showed only shallow copy-number changes and was not associated with
overall survival. Genome-wide, 163 hypermethylated-downregulated and 42 hypomethylated-upregulated genes showed a
significant negative methylation–expression correlation.

**Conclusion:** The deep learning framework identified RAP1GAP2 as a latent, classification-informative gene in head and
neck cancer that is overlooked by DGE analysis. Biological interpretation suggests a possible role through Rap1, MAPK
signalling and Golgi-mediated secretion that now requires experimental validation.

---

## Validation experiments

To validate the computational findings, **qRT-PCR assays** on **head and neck cancer samples** are under way.

- **Collaboration:** a joint project between **BMU (Bangladesh Medical University)** and **SUST (Shahjalal University of Science and Technology)**
- **Ethical clearance:** approved by the BMU and SUST institutional review boards
- **Sample collection:** in progress through the **Department of Otolaryngology – Head & Neck Surgery, BMU**
- **Planned experiments:** qRT-PCR quantification of RAP1GAP2 expression and correlation with clinicopathological features

---

## Repository structure

```
zenodo_manuscript/  Preprint (text source, build script, figures and supplementary tables)
deepprofile/        VAE training and Integrated Gradients (adapted from DeepProfile), result plots,
                    single-gene classifiers and survival analysis (see deepprofile/README.md)
preprocessing/      Dataset download and QC (GSE178537, GSE290057, TCGA-HNSC)
methylation/        TCGA-HNSC promoter methylation × DEG integration (see methylation/README.md)
rap1gap2_analysis/  RAP1GAP2/RAP1GAP multi-omic follow-up: expression (2 cohorts), methylation,
                    CNA, mutations, survival, downstream markers (see rap1gap2_analysis/README.md)
manuscript/         Thesis-era R Markdown and table sources
validation/         BMU–SUST qRT-PCR validation project (see validation/README.md)
archive/            Early test runs, kept for reference
```

---

## Acknowledgements

- The modelling code and the head and neck cancer compendium come from **DeepProfile** (Qiu et al., 2025).
- TCGA data were provided by the TCGA Research Network through the GDC and cBioPortal.

---

## License

Academic / non-commercial use only. Please cite if you use this work.
