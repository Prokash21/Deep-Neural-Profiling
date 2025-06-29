

# Install required packages (if not already installed)
install.packages("flextable")
install.packages("officer")
install.packages("magrittr")

# Load libraries
library(flextable)
library(officer)
library(magrittr)

# Updated R data (18 datasets)
geo_data <- data.frame(
  GEO_Accession = c(
    "GSE37991", "GSE23558", "GSE25099", "GSE10121", "GSE31853",
    "GSE131182", "GSE145272", "GSE217142", "GSE85195", "GSE168227",
    "GSE84805", "GSE30784", "GSE2280", "GSE3524", "GSE6791", 
    "GSE41442", "GSE37371", "GSE23030", "GSE29000"
  ),
  Samples = c(
    "80 (40 tumor + 40 normal)", "31 (27 tumor + 4 normal)", "79 (57 tumor + 22 normal)",
    "41 (35 tumor + 6 normal)", "11 (8 tumor cell lines + 3 normal)",
    "12 (6 paired tumor + normal)", "10 (5 metastatic + 5 non-metastatic)",
    "6 (primary + recurrent tumors)", "49 (34 OSCC + 15 OPL)", "6 paired tumor-normal samples",
    "6 paired tumor-normal samples", "229 total (167 tumor + others)",
    "32 (27 non-metastatic + 5 metastatic)", "20 (16 tumor + 4 normal)",
    "154 (119 tumor + 35 controls)", "55 (45 tumor + 10 normal)",
    "100 (50 tumor + 50 normal)", "30 metastatic tongue OSCC", "50 (40 tumor + 10 normal)"
  ),
  Platform = c(
    "GPL6883 (Illumina HumanRef‑8)", "GPL6480 (Agilent 44K)", "GPL5175 (Affymetrix Exon ST)",
    "Operon Oligoset 4.0", "GPL96/570 (Affymetrix)", "GPL20301 (Illumina HiSeq)",
    "HiSeq 2500 RNA‑seq", "NovaSeq 6000 RNA‑seq", "GPL6480 (Agilent 44K)",
    "Agilent lncRNA microarray", "Agilent lncRNA array", "GPL570 (Affymetrix U133 Plus 2.0)",
    "GPL96 (Affymetrix U133A)", "GPL96 (Affymetrix U133A)", "Affymetrix U133 Plus 2.0",
    "GPL570 (Affymetrix)", "GPL96 (Affymetrix)", "GPL5175 (Affymetrix Exon ST)", "GPL570 (Affymetrix)"
  ),
  Study_Type = c(
    "Expression profiling by array", "Expression profiling by array", "Expression profiling by array",
    "Expression profiling by array", "Expression profiling by array", "Expression profiling by RNA‑seq",
    "Expression profiling by RNA‑seq", "Expression profiling by RNA‑seq", "Expression profiling by array",
    "Expression profiling by array", "Expression profiling by array", "Expression profiling by array",
    "Expression profiling by array", "Expression profiling by array", "Expression profiling by array",
    "Expression profiling by array", "Expression profiling by array", "Expression profiling by array", "Expression profiling by array"
  )
)


# Create the styled flextable
geo_table <- flextable(geo_data) %>%
  theme_vanilla() %>%
  fontsize(size = 10, part = "all") %>%
  autofit() %>%
  set_table_properties(width = 1, layout = "autofit") %>%
  bold(part = "header") %>%
  align(align = "center", part = "all") %>%
  padding(padding = 5)

# Export to Word file
doc <- read_docx()
doc <- body_add_par(doc, "Expression Profiling GEO Datasets for OSCC", style = "heading 1")
doc <- body_add_flextable(doc, value = geo_table)
print(doc, target = "OSCC_GEO_Datasets_Table.docx")

