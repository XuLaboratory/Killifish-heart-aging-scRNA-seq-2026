# Killifish-heart-aging-scRNA-seq-2026
scRNA-seq profiles of hearts from young and old GRZ strain ATK 

## Overview

This repository contains R scripts for analyzing single-cell RNA-seq data from African Turquoise Killifish (*Nothobranchius furzeri*) hearts, comparing young adult (Week 8, W8) versus aging (Week 16, W16) timepoints.

The analysis pipeline consists of 7 steps, covering data integration, annotation, cell proportion analysis, differential expression, trajectory analysis, pathway enrichment, and cell-cell communication.

---

## Scripts & Analysis Pipeline

### `1. Data Integration.R`

- **Purpose**: To integrate week 8 and week 16 heart datasets

- **Key Functions**: 
  - Import week 8 and week 16 gene expression profiles, cell metadata, and matrix files
  - Create Seurat objects for both week 8 and week 16 with quality control
  - Normalize, merge, scale, integrate using Harmony, and cluster the week 8 and week 16 Seurat objects
  
- **Output**:
  - `kf_RIMD.rds`: Integrated Seurat object

---

### `2. Data Annotation.R`

- **Purpose**: To annotate cell clusters for downstream analyses
  
- **Key Functions**:
  - Search the top 50 markers of each cluster
  - Annotate clusters based on known markers curated from the literature
  - Remove clusters containing heterogeneous and ambiguous gene expression profiles
    
- **Outputs**:
  - `kf_simplified.rds`: Annotated Seurat object
  - `stacked_vlnplot.pdf`: Stacked Violin plot visualizing known markers
  - `avg_exp_heatmap.png` : Heatmap plot visualizing known markers
    
---

### `3. Cell proportion analysis.R`

- **Purpose**: To analyze cell type composition changes during cardiac aging (W8 vs W16).

- **Key Functions**:
  - Calculate cell proportions for each cluster in W8 vs W16
  - Create individual comparison plots for each cell type (with 5% Y-axis intervals)
  - Produce multi-panel layouts sorted by total proportion

- **Outputs**:
  - `Individual_*.png`: Per-cluster comparison plots
  - `Combined_5col.png`: 5-column layout
  - `Cell_Proportions.xlsx`: Statistical data table

---

### `4. Differential expression analysis.R`

- **Purpose**: To identify differentially expressed genes between W8 and W16 for each cell type.

- **Key Functions**:
  - Perform Wilcoxon rank-sum test for M, B, EC, VCM, Epi, and T cell cluster
  - Generate volcano plots with labeled top genes
  - Create publication-ready multi-panel volcano plots

- **Outputs**:
  - DEGs_Top50.xlsx: Top 50 DEGs per cluster
  - Volcano_*.pdf/png: Individual volcano plots
  - Merged_DEGs.xlsx: DEGs for merged groups
  - Volcano_Horizontal/Grid.pdf/png: Multi-panel plots

---

### `5. Trajectory analysis.R`

- **Purpose**: To identify cellular transition within VCM (ventricular cardiomyocyte) cluster

- **Key Functions**:
  - Perform the monocle3 package to identify cellular transition within VCM clusters 
  - Generate igraph plots to visualize cellular transition  

- **Outputs**:
  - `plot_list$W8`: Trajectory plot for W8 hearts
  - `plot_list$W16`: Trajectory plot for W16 hearts
 

---

## Requirements

- **R Version**: ≥ 4.3.0

- **Input Data**: `kf_simplified.rds`: Annotated Seurat object containing:
  - 10 cell clusters (M, B, EC, VCM, EC, Epi, T, EPDC, Gran, SMC, Neu)
  - Metadata: `condition` (W8/W16)

---

## Contributors

- Scripts 1, 2, 5: Baul Yoon
- Scripts 3, 4: Edward (Yifeng) Xu

