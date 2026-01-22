# Killifish-heart-aging-scRNA-seq-2026
scRNA-seq profiles of hearts from young and old GRZ strain ATK 

## Overview

This repository contains R scripts for analyzing single-cell RNA-seq data from African Turquoise Killifish (*Nothobranchius furzeri*) hearts, comparing young adult (Week 8, W8) versus aging (Week 16, W16) timepoints.

The analysis pipeline consists of 6 steps, covering data integration, annotation, cell proportion analysis, differential expression, pathway enrichment, and cell-cell communication.

---

## Scripts & Analysis Pipeline

### `1. Data Integration.R`
*(To be completed by colleague)*
- **Purpose**: [Pending description]
- **Key Functions**:
  - ...
- **Outputs**:
  - ...

### `2. Data Annotation.R`
*(To be completed by colleague)*
- **Purpose**: [Pending description]
- **Key Functions**:
  - ...
- **Outputs**:
  - ...

---

### `3. Cell proportion analysis.R`

**Purpose**: To analyze cell type composition changes during cardiac aging (W8 vs W16).

**Key Functions**:
- Calculate cell proportions for each cluster in W8 vs W16
- Create individual comparison plots for each cell type (with 5% Y-axis intervals)
- Produce multi-panel layouts sorted by total proportion

**Outputs**:
- `Individual_*.png`: Per-cluster bar plots (N = number of clusters)
- `Combined_5col.png`: 5-column layout
- `Cell_Proportions.xlsx`: Statistical data table

---

### `4. Differential expression analysis.R`

**Purpose**: To identify differentially expressed genes (DEGs) between W8 and W16 for each cell type.

**Key Functions**:
- Perform Wilcoxon rank-sum test for each cluster
- Generate volcano plots with labeled top genes
- Analyze merged cell type groups (M, B, EC, VCM, Epi, T) for increased statistical power
- Create publication-ready multi-panel volcano plots

**Outputs**:
- `DEGs_Top50.xlsx`: Top 50 DEGs per cluster (ranked by |log2FC|)
- `Volcano_*.pdf/png`: Individual volcano plots for each cluster
- `Merged_DEGs.xlsx`: DEGs for merged groups
- `Volcano_Merged_*.pdf/png`: Volcano plots for merged groups
- `Volcano_Grid.pdf/png`: 2×3 layout for publication
- `All_DEGs.xlsx`: Complete DEG results for all clusters

---

### `5. Pathway enrichment gsea.R`

**Purpose**: To perform Gene Set Enrichment Analysis (GSEA) using human orthologs to identify pathway-level changes.

**Key Functions**:
- Convert killifish gene names to human orthologs
- Map gene symbols to Entrez IDs for pathway analysis
- Run gseGO (Gene Ontology Biological Process) for major cell type groups
- Generate ridge plots showing enrichment score distributions
- Extract core enrichment genes driving pathway significance

**Outputs**:
- `GSEA_Results.xlsx`: Complete GSEA results (GO terms, NES, p-values)
- `GSEA_RidgePlots.pdf`: Ridge plots for all groups
- `GO_Gene_Lists.xlsx`: Top 10 core enrichment genes per pathway
- `RidgePlot_Gene_Summary.xlsx`: Top 5 genes per pathway (for ridge plot interpretation)

---

### `6. Cell-cell communication.R`
*(To be completed by colleague)*
- **Purpose**: [Pending description]
- **Key Functions**:
  - ...
- **Outputs**:
  - ...

---

## Requirements

**R Version**: ≥ 4.3.0

**Input Data**:
- `kf_simplified.rds`: Preprocessed Seurat object containing:
  - 10 cell clusters (M, B, EC, VCM, EC, Epi, T, EPDC, Gran, SMC, Neu)
  - Metadata: `condition` (W8/W16)

---
**Contributors**:
- Scripts 1, 2, 6: Baul Yoon
- Scripts 3, 4, 5: Edward (Yifeng) Xu

