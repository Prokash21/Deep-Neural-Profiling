# Zenodo preprint package

Preprint: *Deep Neural Profiling Reveals RAP1GAP2 as a Latent Regulator in Head and Neck Cancer: Integrative Transcriptomic and Multi-omic Analysis*.

Based on the undergraduate thesis in `manuscript/`, with the TCGA-HNSC methylation and RAP1GAP2 multi-omic analyses added.

## Files

| File | Content |
|---|---|
| `manuscript.pdf` | Manuscript for upload (23 pages; 9 figures, 3 tables, Supplementary Figure S1) |
| `manuscript.docx` | Editable Word version (same content) |
| `supplementary/Supplementary_Tables_S1-S10.xlsx` | Supplementary Tables S1–S10 |
| `figures/` | Figures 1–9 and S1 (PNG, 300 dpi) |
| `manuscript_content.py` | Text of the manuscript (single source for the PDF and DOCX) |
| `build_manuscript.py` | Rebuilds `manuscript.docx` and `manuscript.html` from `manuscript_content.py` |

To rebuild after editing the text:

```bash
pip install python-docx
python build_manuscript.py
"/Applications/Google Chrome.app/Contents/MacOS/Google Chrome" --headless --no-pdf-header-footer \
    --virtual-time-budget=10000 --print-to-pdf="$PWD/manuscript.pdf" "file://$PWD/manuscript.html"
```

## Sources of every result

| Manuscript item | Source in this repository |
|---|---|
| Figures 1–7, 9; Tables 2–3 | `deepprofile/illustrator_figures/` and the thesis (`manuscript/`) |
| VAE / IG / classifier settings | `deepprofile/MODEL_TRAININGS/`, `deepprofile/Single_Gene_Model/ML/` |
| Figure 8A, 8C, 8D; Tables S6–S9 | `rap1gap2_analysis/` |
| Figure 8B, Figure S1; Tables S2–S5, S10 | `methylation/results/` |
| Table 1, Table S1 | NCBI GEO metadata (`manuscript/Table/Table_2.1_GEO_datasets.csv`) |

## Suggested Zenodo metadata

- **Resource type:** Publication → Preprint
- **Title:** as above
- **Creators:** Debnath, Joy Prokash; Hossain, Tanvir; Rahman, Papia (all: Department of Biochemistry and Molecular Biology, Shahjalal University of Science and Technology)
- **Description:** the abstract from the manuscript
- **Keywords:** head and neck cancer; variational autoencoder; DeepProfile; Integrated Gradients; RAP1GAP2; DNA methylation; TCGA-HNSC; multi-omics
- **License:** Creative Commons Attribution 4.0 International (CC BY 4.0)
- **Related identifiers:**
  - *is supplemented by* https://github.com/Prokash21/Deep-Neural-Profiling (software)
  - *references* 10.1038/s41551-024-01290-8 (DeepProfile, Qiu et al. 2025)
- **Files to upload:** `manuscript.pdf`, `supplementary/Supplementary_Tables_S1-S10.xlsx` (optionally the `figures/` PNGs)
