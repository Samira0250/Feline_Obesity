# Feline Adipose Tissue Transcriptomics

**Species:** *Felis catus* (domestic cat)  
**Tissue:** Abdominal and subcutaneous adipose  
**Design:** Paired repeated measures — 8 cats × 3 timepoints  
**Analyses:** Differential gene expression · GO enrichment · Immune cell deconvolution

---

## Overview

This repository contains all bioinformatics code used to analyze bulk RNA-seq data from feline adipose tissue collected across three dietary timepoints: **Baseline**, **Developing Obesity** (Middle), and **Prolonged Obesity** (End). Two adipose depots were profiled — abdominal and subcutaneous — allowing both temporal and anatomical comparisons.

The pipeline covers two independent analyses:

1. **Differential Gene Expression (DGE)** — identifying genes that change across timepoints using DESeq2, visualized with volcano plots, heatmaps, and GO Biological Process enrichment.
2. **Immune Cell Deconvolution** — estimating immune cell-type proportions in bulk RNA-seq data using a human single-cell reference and NNLS regression, validated against CIBERSORTx.

> Associated manuscript: *Chronic Obesity in Domestic Cats Does Not Induce a Classic Inflammatory Response*  
> GEO accession: *PRJNA940803*

---

## Repository Structure

```
feline-adipose-rnaseq/
│
├── dge/                          # Differential gene expression analysis
│   ├── scripts/
│   │   ├── 01_dge_volcanos_heatmaps_GO.R   # MID vs BSL, END vs MID
│   │   └── 02_end_vs_bsl.R                 # END vs BSL
│   ├── data/                     # Input files (not tracked — see below)
│   ├── plots2/                   # All DGE outputs
│   ├── .gitignore
│   └── README.md
│
├── deconvolution/                # Immune cell deconvolution analysis
│   ├── scripts/
│   │   └── adipose_immune_deconvolution.R
│   ├── data/                     # Input files (not tracked — see below)
│   ├── deconv_output/            # All deconvolution outputs
│   ├── .gitignore
│   └── README.md
│
└── README.md                     # This file
```

Each subfolder has its own `README.md` with full installation, usage, and output documentation.

---

## Experimental Design

| Factor | Levels |
|---|---|
| Depot | Abdominal, Subcutaneous |
| Timepoint | Baseline, Developing Obesity (Middle), Prolonged Obesity (End) |
| Animals | 8 cats (paired, repeated measures) |
| DESeq2 model | `~ Animal + Timepoint` |
| Comparisons | MID vs BSL · END vs MID · END vs BSL |

---

## Analysis 1 — Differential Gene Expression

**Scripts:** `dge/scripts/01_dge_volcanos_heatmaps_GO.R` and `02_end_vs_bsl.R`

DESeq2 results are read from pre-computed output sheets. Log2 fold changes are flipped so that positive values represent upregulation at the *later* timepoint. Genes are classified as differentially expressed at **|log2FC| > 1** and **padj < 0.05**.

**Key outputs:**

| Output | Description |
|---|---|
| `volcano_2x2.pdf` | Combined 2×2 volcano grid (AB + SQ, MID vs BSL + END vs MID) |
| `volcano_END_vs_BSL.pdf` | Side-by-side volcanos for END vs BSL |
| `heatmaps.pdf` | Top 50 DEGs per comparison, Z-scored pseudo counts |
| `go_dotplots.pdf` | GO Biological Process enrichment (UP, DOWN, combined) |
| `DEG_up_down_tables.xlsx` | Filtered DEG tables per comparison |
| `GO_results.xlsx` | GO enrichment results per comparison |

**R packages:** `DESeq2` · `ggplot2` · `ggrepel` · `pheatmap` · `clusterProfiler` · `org.Hs.eg.db` · `patchwork`

> Note: GO enrichment uses human gene symbols (`org.Hs.eg.db`) as a proxy for feline annotation, which is standard practice given the limited functional annotation of the *Felis catus* genome.

---

## Analysis 2 — Immune Cell Deconvolution

**Script:** `deconvolution/scripts/adipose_immune_deconvolution.R`

Bulk RNA-seq counts are normalized to TPM, then cat Ensembl IDs are mapped to human gene symbols via case-insensitive string matching. Cell-type proportions are estimated by NNLS regression against a human immune single-cell reference — algorithmically equivalent to CIBERSORTx. Output files are formatted for direct upload to the CIBERSORTx portal for independent validation.

**Cell types detected (12):**

| Category | Cell Types |
|---|---|
| Myeloid — Tissue Resident | Anti-inflammatory Macrophage (Mac1), Proinflammatory Macrophage (Mac2+3) |
| Myeloid — Circulating | Monocyte (Mono2) |
| Myeloid — Antigen Presenting | Dendritic cells (cDC1+2), Plasmacytoid DC, ASDC |
| Myeloid — Effector | Mast cells, Neutrophil |
| Lymphoid | NK cells, Plasma cell, Th1 cell, Th2 cell |

**Key outputs:**

| Output | Description |
|---|---|
| `adipose_immune_deconvolution_results.csv` | Per-sample proportions + metadata |
| `immune_summary_proportion_table.csv` | Mean % per cell type, depot, and timepoint |
| `immune_statistical_tests_*.csv` | Paired Wilcoxon tests with FDR correction |
| `IMMUNE_publication_*.pdf/.png` | Publication-ready figures (3 panels) |
| `CIBERSORTx_*.txt` | Files for CIBERSORTx web validation |

**R packages:** `Seurat` · `nnls` · `ggplot2` · `pheatmap` · `viridis` · `clusterProfiler` · `DESeq2`

---

## Requirements

**R version ≥ 4.2 | Bioconductor ≥ 3.16**

```r
# CRAN
install.packages(c("Seurat", "nnls", "openxlsx", "tidyverse", "ggrepel",
                   "patchwork", "pheatmap", "RColorBrewer", "reshape2", "viridis"))

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "biomaRt", "SingleCellExperiment",
                       "clusterProfiler", "org.Hs.eg.db"))
```

---

## Data Availability

Raw data files are not tracked in this repository due to size. Input files should be placed in the `data/` folder of each analysis directory before running scripts. See individual `README.md` files for the complete list of required files.

Processed data and raw counts are available at:  
**GEO:** *[PRJNA940803]*

---

## Citation

If you use this code, please cite:

> *[DOI to paper when published]*

---

## Author

**Samira Salimiyan**  
Auburn University  
Contact: [szs0394@auburn.edu]
