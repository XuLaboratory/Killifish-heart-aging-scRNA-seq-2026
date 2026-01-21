# Killifish-heart-aging-scRNA-seq-2026
scRNA-seq profiles of hearts from young and old GRZ strain ATK 

## Overview

This repository contains R scripts for analyzing single-cell RNA-seq data from African Turquoise Killifish (*Nothobranchius furzeri*) hearts, comparing young adult (Week 8, W8) versus aging (Week 16, W16) timepoints.

The analysis pipeline consists of 6 steps, covering data integration, annotation, cell proportion analysis, differential expression, pathway enrichment, and cell-cell communication.

---

## Scripts & Analysis Pipeline

### 1. `1. Data Integration.R`
*(To be completed by colleague)*
- **Purpose**: [Pending description]
- **Key Functions**:
  - ...
- **Outputs**:
  - ...

### 2. `2. Data Annotation.R`
*(To be completed by colleague)*
- **Purpose**: [Pending description]
- **Key Functions**:
  - ...
- **Outputs**:
  - ...

---

### 3. `3. Cell proportion analysis.R`

**Purpose**: Analyzes cell type composition changes during cardiac aging (W8 vs W16).

**Key Functions**:
- Calculates cell proportions for each cluster in W8 vs W16
- Creates individual comparison plots for each cell type (with 5% Y-axis intervals)
- Produces multi-panel layouts sorted by total proportion

**Outputs**:
- `Q1_Individual_*.png`: Per-cluster bar plots (N = number of clusters)
- `Q1_Combined_5col.png`: 5-column layout
- `Q1_Combined_4col.png`: 4-column layout
- `Q1_Cell_Proportions.xlsx`: Statistical data table

**Runtime**: ~2-5 minutes

---

### 4. `4. Differential expression analysis.R`

**Purpose**: Identifies differentially expressed genes (DEGs) between W8 and W16 for each cell type.

**Key Functions**:
- Performs Wilcoxon rank-sum test for each cluster
- Generates volcano plots with labeled top genes
- Analyzes merged cell type groups (M, B, EC, VCM, Epi, T) for increased statistical power
- Creates publication-ready multi-panel volcano plots

**Outputs**:
- `Q2_DEGs_Top50.xlsx`: Top 50 DEGs per cluster (ranked by |log2FC|)
- `Q2_Volcano_*.pdf/png`: Individual volcano plots for each cluster
- `Q2_Merged_DEGs.xlsx`: DEGs for merged groups
- `Q2_Volcano_Merged_*.pdf/png`: Volcano plots for merged groups
- `Q2_Volcano_Horizontal.pdf/png`: 1×6 layout for publication
- `Q2_Volcano_Grid.pdf/png`: 2×3 layout for publication
- `Q2_All_DEGs.xlsx`: Complete DEG results for all clusters

**Runtime**: ~10-30 minutes (depending on cell counts)

---

### 5. `5. Pathway enrichment gsea.R`

**Purpose**: Performs Gene Set Enrichment Analysis (GSEA) using human orthologs to identify pathway-level changes.

**Key Functions**:
- Converts killifish gene names to human orthologs
- Maps gene symbols to Entrez IDs for pathway analysis
- Runs gseGO (Gene Ontology Biological Process) for major cell type groups
- Generates ridge plots showing enrichment score distributions
- Extracts core enrichment genes driving pathway significance

**Outputs**:
- `Q3_GSEA_Results.xlsx`: Complete GSEA results (GO terms, NES, p-values)
- `Q3_GSEA_RidgePlots.pdf`: Ridge plots for all groups
- `Q3_GO_Gene_Lists.xlsx`: Top 10 core enrichment genes per pathway
- `Q3_RidgePlot_Gene_Summary.xlsx`: Top 5 genes per pathway (for ridge plot interpretation)

**Runtime**: ~5-15 minutes

---

### 6. `6. Cell-cell communication.R`
*(To be completed by colleague)*
- **Purpose**: [Pending description]
- **Key Functions**:
  - ...
- **Outputs**:
  - ...

---

## Requirements

**R Version**: ≥ 4.3.0

**CRAN Packages**:
- Seurat, dplyr, ggplot2, writexl, readxl, ggrepel, patchwork, RColorBrewer, tidyr

**Bioconductor Packages**:
- clusterProfiler, org.Hs.eg.db, enrichplot

**Input Data**:
- `kf_simplified.rds`: Preprocessed Seurat object containing:
  - 10 merged cell types (VCM, EC, Epi, EPDC, M, B, T, Gran, SMC, Neu)
  - RBC cells already removed
  - Metadata: `condition` (W8/W16)

---

## Usage

**Note**: Scripts 3, 4, and 5 are standalone and include their own environment setup.

**Steps**:
1. Open the desired script and modify `base_dir` to point to your data directory
2. Run scripts sequentially:
```r
# --- Phase 1: Integration & Annotation ---
# source("1. Data Integration.R")
# source("2. Data Annotation.R")

# --- Phase 2: Analysis (Edward) ---
source("3. Cell proportion analysis.R")
source("4. Differential expression analysis.R")
source("5. Pathway enrichment gsea.R")

# --- Phase 3: Communication ---
# source("6. Cell-cell communication.R")
```

---

## Output Structure
```
results/
├── Q1_*.png/xlsx              # Cell proportion results (Script 3)
├── Q2_*.pdf/png/xlsx          # DEG and volcano plot results (Script 4)
└── Q3_*.pdf/xlsx              # GSEA pathway results (Script 5)
```

---

## Citation & Contact

**Citations**:
- Seurat: Hao et al. (2021) *Cell*
- clusterProfiler: Yu et al. (2012) *OMICS*

**Contributors**:
- Scripts 3, 4, 5: Edward (Yifeng Xu) - Bioinformatics Volunteer, Maryland
- Scripts 1, 2, 6: [Colleague Name/Info]
