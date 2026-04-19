# 🧬 Parkinson's Disease — RNA-Seq Exploratory Analysis & Biomarker Discovery

[![R](https://img.shields.io/badge/R-4.3.3-276DC3?style=flat&logo=r&logoColor=white)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-ComplexHeatmap-87CEEB)](https://bioconductor.org/packages/ComplexHeatmap/)
[![Course](https://img.shields.io/badge/Course-Multiomics%20Data%20Analysis-brightgreen)](.)

> Identification of the most variable genes as potential biomarkers for Parkinson's Disease using microarray gene expression data, as part of the **Multiomics Data Analysis** course — Bioinformatics Diploma Program.

---

## 📋 Table of Contents

- [Background](#-background)
- [Dataset](#-dataset)
- [Pipeline Overview](#-pipeline-overview)
- [Results Summary](#-results-summary)
- [Repository Structure](#-repository-structure)
- [Requirements](#-requirements)
- [How to Run](#-how-to-run)
- [Output Files](#-output-files)
- [Key Findings](#-key-findings)
- [References](#-references)

---

## 🧠 Background

**Parkinson's Disease (PD)** is the second most common neurodegenerative disorder, characterized by progressive loss of dopaminergic neurons in the substantia nigra. Early diagnosis remains a major clinical challenge.

This project applies **exploratory RNA-Seq analysis** to:
1. Assess data quality across all samples
2. Identify the most transcriptionally variable genes between PD patients and healthy controls
3. Visualize expression patterns using dimensionality reduction (PCA) and clustering (heatmaps)

---

## 📊 Dataset

| Property | Value |
|---|---|
| Platform | Microarray (Affymetrix) |
| Samples | 19 total (10 Ctrl, 9 PD) |
| Features | 22,283 probe-gene pairs |
| Expression scale | log2 (pre-normalized) |
| Age range | 68–89 years |

### Sample breakdown

| Group | N | Gender |
|---|---|---|
| Control (Ctrl) | 10 | 3F / 7M |
| Parkinson's Disease (PD) | 9 | 3F / 6M |

---

## 🔬 Pipeline Overview

```
Raw expression matrix (22,283 × 19)
          │
          ▼
  ┌───────────────────┐
  │  1. Load & QC     │  → Validate sample names, dimensions
  └────────┬──────────┘
           │
           ▼
  ┌───────────────────┐
  │  2. EDA           │  → Histogram · Density plot per sample
  └────────┬──────────┘
           │
           ▼
  ┌───────────────────┐
  │  3. PCA           │  → 2D (ggplot2) · 3D (plotly)
  └────────┬──────────┘
           │
           ▼
  ┌───────────────────┐
  │  4. Row Variance  │  → rowVars() → top 100 genes
  └────────┬──────────┘
           │
           ▼
  ┌──────────────────────────────────┐
  │  5. Heatmaps                     │
  │   • Absolute expression          │
  │   • Z-score (relative)  [BONUS]  │
  └──────────────────────────────────┘
           │
           ▼
  Top 100 variable genes → Candidate PD biomarkers
```

---

## 📈 Results Summary

### Quality Control
- All 19 samples passed QC — density curves are tightly overlapping
- Expression values follow an expected right-skewed distribution (log2 scale: 6–14)

### PCA
| Component | Variance Explained |
|---|---|
| PC1 | 25.4% |
| PC2 | 17.1% |
| PC3 | 11.8% |
| **Total (PC1–3)** | **54.3%** |

- Partial separation between PD and Ctrl groups on PC1 → disease-driven signal present
- No clear gender-based clustering → disease effect dominates

### Top 10 Most Variable Genes

| Rank | Probe | Gene | Variance | Biological Role |
|---|---|---|---|---|
| 1 | 204141_at | **TUBB2A** | 2.147 | Tubulin — cytoskeletal integrity |
| 2 | 203282_at | **GBE1** | 2.037 | Glycogen metabolism |
| 3 | 200799_at | — | 1.991 | — |
| 4 | 202581_at | — | 1.848 | — |
| 5 | 200633_at | **UBB** | 1.824 | Ubiquitin — protein degradation |
| 6 | 202482_x_at | **RANBP1** | 1.732 | Nuclear transport |
| 7 | 200863_s_at | **RAB11A** | 1.650 | Vesicle trafficking |
| 8 | 209118_s_at | **TUBA1A** | 1.643 | Tubulin alpha — axonal transport |
| 9 | 208609_s_at | — | 1.638 | — |
| 10 | 208845_at | **VDAC3** | 1.613 | Mitochondrial membrane channel |

---

## 🗂️ Repository Structure

```
parkinson-rnaseq-analysis/
│
├── data/
│   ├── Parkinson_exp.txt             # Expression matrix (22,283 × 19)
│   └── Parkinson_phenotable.txt      # Sample metadata
│
├── scripts/
│   └── 01_parkinson_analysis.R       # Main analysis script
│
├── results/
│   ├── top100_variable_genes.csv     # Top 100 genes ranked by variance
│   └── plots/
│       ├── 01_Histogram_Expression.pdf
│       ├── 02_Density_Plot_per_Sample.pdf
│       ├── 03_PCA_2D.pdf
│       ├── 04_PCA_3D.html            # Interactive 3D plot
│       ├── 05_Heatmap_Absolute_Expression.pdf
│       └── 06_Heatmap_Zscore.pdf
│
├── docs/
│   └── graphical_abstract.svg        # Workflow diagram
│
├── .gitignore
├── LICENSE
└── README.md
```

---

## ⚙️ Requirements

### R version
```
R >= 4.3.0
```

### CRAN packages
```r
install.packages(c("ggplot2", "plotly", "matrixStats", "htmlwidgets", "circlize"))
```

### Bioconductor packages
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
```

---

## 🚀 How to Run

### 1. Clone the repository
```bash
git clone https://github.com/<your-username>/parkinson-rnaseq-analysis.git
cd parkinson-rnaseq-analysis
```

### 2. Add your data files
Place the two data files inside the `data/` folder:
```
data/Parkinson_exp.txt
data/Parkinson_phenotable.txt
```

### 3. Run the analysis
Open RStudio, set working directory to `scripts/`, then run:
```r
setwd("scripts/")
source("01_parkinson_analysis.R")
```

Or from terminal:
```bash
Rscript scripts/01_parkinson_analysis.R
```

All output PDFs and CSVs will be saved automatically to `results/`.

---

## 📁 Output Files

| File | Description |
|---|---|
| `01_Histogram_Expression.pdf` | Distribution of all expression values |
| `02_Density_Plot_per_Sample.pdf` | Per-sample expression density curves |
| `03_PCA_2D.pdf` | PCA scatter plot with confidence ellipses |
| `04_PCA_3D.html` | Interactive 3D PCA (open in browser) |
| `05_Heatmap_Absolute_Expression.pdf` | Heatmap of top 100 genes (raw values) |
| `06_Heatmap_Zscore.pdf` | Heatmap of top 100 genes (z-score) |
| `top100_variable_genes.csv` | Ranked table of top 100 variable genes |

---

## 🔍 Key Findings

The analysis revealed a set of highly variable genes between PD and healthy controls, with strong biological relevance to PD pathology:

- **TUBB2A & TUBA1A** — Tubulin genes involved in cytoskeletal dynamics and axonal transport, processes known to be disrupted in PD
- **UBB (Ubiquitin B)** — Core component of the ubiquitin-proteasome system; implicated in Lewy body formation — a hallmark of PD
- **VDAC3** — Mitochondrial outer membrane channel; mitochondrial dysfunction is a central mechanism in PD neurodegeneration
- **RAB11A** — Involved in vesicle recycling and dopamine receptor trafficking

These genes represent strong candidates for further investigation as diagnostic biomarkers or therapeutic targets.

---

## 📚 References

1. Hamed M. et al. (2018). *A workflow for the integrative transcriptomic description of molecular pathology and the suggestion of normalizing compounds, exemplified by Parkinson's disease.* Scientific Reports, 8, 7937. https://doi.org/10.1038/s41598-018-25754-5

2. Gu Z. et al. (2016). *Complex heatmaps reveal patterns and correlations in multidimensional genomic data.* Bioinformatics, 32(18), 2847–2849.

3. R Core Team (2024). *R: A Language and Environment for Statistical Computing.* https://www.r-project.org/

---

## 👤 Author

**Ibrahim (Bembo)**
Bioinformatician | Field Application Specialist
MSc Bioinformatics — University of Sadat City
Bioinformatics Diploma — Multiomics Data Analysis Track

---

## 📄 License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.
