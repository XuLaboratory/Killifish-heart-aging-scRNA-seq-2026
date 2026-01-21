library(CellChat)
library(patchwork)
library(Seurat)
library(biomaRt)
library(Matrix)

kf = readRDS("kf_simplified.rds") 

# --- Mouse_Mart generation and ATK homology gene search ---
CellChatDB <- CellChatDB.mouse
kf$cluster <- factor(Idents(kf))

kf_list <- SplitObject(kf, split.by = "condition")

kf_W8  <- kf_list$W8
kf_W16 <- kf_list$W16

expr <- GetAssayData(kf, assay = "RNA", layer = "data")
gene <- rownames(expr)[
  grepl("^[a-zA-Z]", rownames(expr)) &        
    !grepl("^KEG", rownames(expr)) &             
    !grepl("^trna", rownames(expr), ignore.case = TRUE)
]

mouse_mart <- useEnsembl(
  biomart = "genes",
  dataset = "mmusculus_gene_ensembl",
  mirror  = "useast"
)

homology_map <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "nfurzeri_homolog_associated_gene_name"
  ),
  mart = mouse_mart
)

gene_map <- getBM(
  attributes = c(
    "ensembl_gene_id", "mgi_symbol"),
  mart = mouse_mart
)

homology_map <- getBM(attributes = c("ensembl_gene_id", "nfurzeri_homolog_associated_gene_name"), mart = mouse_mart)
gene_map <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"), mart = mouse_mart)

mapping_mouse <- homology_map |>
  inner_join(gene_map, by = "ensembl_gene_id") |>
  filter(nfurzeri_homolog_associated_gene_name %in% gene, mgi_symbol != "")

gene.map <- mapping_mouse$mgi_symbol
names(gene.map) <- mapping_mouse$nfurzeri_homolog_associated_gene_name

new.genes <- gene.map[rownames(kf)]

keep <- !is.na(new.genes) & new.genes != ""

kf_mouse <- kf[keep, ]
rownames(kf_mouse) <- new.genes[keep]

expr <- GetAssayData(kf_mouse, assay = "RNA", layer = "data")
genes <- rownames(expr)

expr <- GetAssayData(
  kf_mouse,
  assay = "RNA",
  layer = "data"   # log-normalized, correct for CellChat
)

genes <- rownames(expr)
gene_levels <- unique(genes)
gene_index <- match(genes, gene_levels)

collapse_mat <- sparseMatrix(
  i = gene_index,
  j = seq_along(gene_index),
  x = 1,
  dims = c(length(gene_levels), length(gene_index))
)

expr_collapsed <- collapse_mat %*% expr

rownames(expr_collapsed) <- gene_levels
kf_mouse_collapsed <- CreateSeuratObject(
  counts = expr_collapsed,
  meta.data = kf_mouse@meta.data
)
kf_mouse_collapsed <- SetAssayData(
  object   = kf_mouse_collapsed,
  assay    = "RNA",
  layer    = "data",
  new.data = expr_collapsed
)

Idents(kf_mouse_collapsed) <- kf_mouse$cluster

mouse.lr.genes <- unique(c(CellChatDB.mouse$interaction$ligand, CellChatDB.mouse$interaction$receptor))
expr.mouse <- GetAssayData(kf_mouse_collapsed, layer = "data")

meta <- data.frame(
  celltype = kf_mouse_collapsed$cluster,
  row.names = colnames(kf_mouse_collapsed)
)

kf_W8  <- subset(kf, subset = condition == "W8")
kf_W16 <- subset(kf, subset = condition == "W16")

map_and_collapse <- function(seurat_obj, gene_map) {
  expr <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")
  
  new.genes <- gene_map[rownames(expr)]
  keep <- !is.na(new.genes) & new.genes != ""
  expr <- expr[keep, ]
  genes <- new.genes[keep]
  
  gene_levels <- unique(genes)
  gene_index  <- match(genes, gene_levels)
  
  collapse_mat <- sparseMatrix(
    i = gene_index,
    j = seq_along(gene_index),
    x = 1,
    dims = c(length(gene_levels), length(gene_index))
  )
  
  expr_collapsed <- collapse_mat %*% expr
  rownames(expr_collapsed) <- gene_levels
  
  out <- CreateSeuratObject(
    counts = expr_collapsed,
    meta.data = seurat_obj@meta.data
  )
  
  out <- SetAssayData(
    object   = out,
    assay    = "RNA",
    layer    = "data",
    new.data = expr_collapsed
  )
  
  Idents(out) <- seurat_obj$cluster
  return(out)
}

# --- W8 Hearts ---

kf_mouse_W8  <- map_and_collapse(kf_W8, gene.map)

expr_W8 <- GetAssayData(kf_mouse_W8, layer = "data")
meta_W8 <- data.frame(celltype = kf_mouse_W8$cluster, row.names = colnames(kf_mouse_W8))
cellchat_W8 <- createCellChat(object = expr_W8, meta = meta_W8, group.by = "celltype")

cellchat_W8@DB <- CellChatDB.mouse
cellchat_W8 <- subsetData(cellchat_W8)
cellchat_W8 <- identifyOverExpressedGenes(cellchat_W8)
cellchat_W8 <- identifyOverExpressedInteractions(cellchat_W8)
cellchat_W8 <- computeCommunProb(cellchat_W8)
cellchat_W8 <- filterCommunication(cellchat_W8, min.cells = 10)
cellchat_W8 <- computeCommunProbPathway(cellchat_W8)
cellchat_W8 <- aggregateNet(cellchat_W8)
cellchat_W8 <- netAnalysis_computeCentrality(cellchat_W8, slot.name = "netP")

netAnalysis_signalingRole_heatmap(cellchat_W8, slot.name = "netP", pattern = "outgoing")
netAnalysis_signalingRole_heatmap(cellchat_W8, slot.name = "netP", pattern = "incoming")

# --- W16 Hearts --- 

kf_mouse_W16 <- map_and_collapse(kf_W16, gene.map)

expr_W16 <- GetAssayData(kf_mouse_W16, layer = "data")
meta_W16 <- data.frame(celltype = kf_mouse_W16$cluster, row.names = colnames(kf_mouse_W16))

cellchat_W16 <- createCellChat(object = expr_W16, meta = meta_W16, group.by = "celltype")
cellchat_W16@DB <- CellChatDB.mouse
cellchat_W16 <- subsetData(cellchat_W16)
cellchat_W16 <- identifyOverExpressedGenes(cellchat_W16)
cellchat_W16 <- identifyOverExpressedInteractions(cellchat_W16)
cellchat_W16 <- computeCommunProb(cellchat_W16)
cellchat_W16 <- filterCommunication(cellchat_W16, min.cells = 10)
cellchat_W16 <- computeCommunProbPathway(cellchat_W16)
cellchat_W16 <- aggregateNet(cellchat_W16)
cellchat_W16 <- netAnalysis_computeCentrality(cellchat_W16)

netAnalysis_signalingRole_heatmap(cellchat_W16, slot.name = "netP", pattern = "outgoing")
netAnalysis_signalingRole_heatmap(cellchat_W16, slot.name = "netP", pattern = "incoming")

# --- Merged --- 
cellchat_merged <- mergeCellChat(list(W8 = cellchat_W8, W16 = cellchat_W16), add.names = c("W8","W16"))
cellchat_list <- list(W8 = cellchat_W8, W16 = cellchat_W16)
cellchat_merged <- mergeCellChat(cellchat_list, add.names = names(cellchat_list))
netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE, measure = "weight")

p <- rankNet(cellchat_merged, mode = "comparison", stacked = TRUE)

p + scale_fill_manual(
  values = c(
    "W8"  = "#377EB8",  
    "W16" = "#E41A1C" 
  )
) +
  scale_fill_manual(
    values = c("W8" = "#377EB8", "W16" = "#E41A1C")
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 14),
    axis.text.y  = element_text(size = 14),
    axis.text.x  = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    plot.title   = element_text(size = 16))
