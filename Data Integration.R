library(Seurat)
library(Matrix)
library(ggplot2)
library(harmony)

# --- Paths to DGE_filtered outputs ---
path_w8   <- "/research/labs/biochem/xux/m310533/seq/files/kf_test/kf_comb/W8/DGE_filtered/" ## Hide /research ~ /kf_comb
path_w16 <- "/research/labs/biochem/xux/m310533/seq/files/kf_test/kf_comb/W16/DGE_filtered/" ## Hide /research ~ /kf_comb

# --- W8 ---
counts_w8 <- readMM(file = paste0(path_w8, "w8_count_matrix.mtx"))
counts_w8 <- t(counts_w8)
dim(counts_w8)
genes_w8  <- read.csv(paste0(path_w8, "w8_all_genes.csv"), header = TRUE)
cells_w8  <- read.csv(paste0(path_w8, "w8_cell_metadata.csv"), header = TRUE)
rownames(counts_w8) = genes_w8$gene_name
colnames(counts_w8) = cells_w8$bc_wells
counts_w8 <- counts_w8[!duplicated(rownames(counts_w8)), ]

w8 <- CreateSeuratObject(counts = counts_w8, project = "W8",
                         min.cells = 3, min.features = 200)
w8$condition <- "W8"

# --- W16 ---
counts_w16 <- readMM(file = paste0(path_w16, "w16_count_matrix.mtx"))
counts_w16 <- t(counts_w16)
dim(counts_w16)
genes_w16  <- read.csv(paste0(path_w16, "w16_all_genes.csv"), header = TRUE)
cells_w16  <- read.csv(paste0(path_w16, "w16_cell_metadata.csv"), header = TRUE)
rownames(counts_w16) = genes_w16$gene_name
colnames(counts_w16) = cells_w16$bc_wells
counts_w16 <- counts_w16[!duplicated(rownames(counts_w16)), ]

w16 <- CreateSeuratObject(counts = counts_w16, project = "W16",
                          min.cells = 3, min.features = 200)
w16$condition <- "W16"

# --- Additional QC filtering ---
w8 <- subset(w8, subset = nFeature_RNA > 100 & nFeature_RNA < 5000)
w16 <- subset(w16, subset = nFeature_RNA > 100 & nFeature_RNA < 5000)

# --- Normalize and find features ---
w8 <- NormalizeData(w8)
w8 <- FindVariableFeatures(w8, selection.method = "vst", nfeatures = 3000)
w16 <- NormalizeData(w16)
w16 <- FindVariableFeatures(w16, selection.method = "vst", nfeatures = 3000)

# --- Merge ---
kf <- merge(w8, y = w16, add.cell.ids = c("W8", "W16"))

# --- Scale, PCA ---
kf <- ScaleData(kf)
kf <- RunPCA(kf, npcs = 50)

# --- Harmony integration ---
kf <- RunHarmony(kf, group.by.vars = "condition", plot_convergence = TRUE, theta = 10)

# --- UMAP / clustering ---
kf <- RunUMAP(kf, reduction = "harmony", dims = 1:25, spread = 2, min.dist = 0.05, n.neighbors = 30)
kf <- FindNeighbors(kf, reduction = "harmony", dims = 1:25)
kf <- FindClusters(kf, resolution = 0.6) 

# --- Visualization ---
DimPlot(kf, reduction = "umap", group.by = "condition", pt.size = 0.5)
saveRDS(kf, "kf_RIMD.rds")


