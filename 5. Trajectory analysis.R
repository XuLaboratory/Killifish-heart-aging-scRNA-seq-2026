# ==============================================================================
# Trajectory analysis - Killifish Heart Aging Study
# 
# Purpose: To identify cellular transition within VCM (ventricular cardiomyocyte) clusters 
#          
# Input: kf_simplified.rds
#
# Output: 
#   - plot_list$W8
#   - plot_list$W16
#   
# Author: Baul Yoon
# ==============================================================================

library(Seurat)
library(monocle3)
library(SeuratWrappers)
library(igraph)
library(ggplot2)

kf = readRDS("kf_simplified.rds") 

conditions <- c("W8", "W16")

get_earliest_principal_node <- function(cds_cm, cluster_name) {
  cell_ids <- which(colData(cds_cm)$cluster == cluster_name)
  closest_vertex <- cds_cm@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  
  tab <- table(closest_vertex[cell_ids])
  best_vertex <- names(tab)[which.max(tab)]  
  
  g <- principal_graph(cds_cm)[["UMAP"]]
  V(g)$name[as.numeric(best_vertex)]
}

cds_list <- list()
plot_list <- list()

for (cond in conditions) {
  
  vcm <- subset(kf, idents = c("VCM1", "VCM2", "VCM3"))
  vcm <- subset(vcm, subset = condition == cond)
  
  cds_cm <- as.cell_data_set(vcm)
  reducedDims(cds_cm)$UMAP <- Embeddings(vcm, "umap")
  
  colData(cds_cm)$cluster <- Idents(vcm)
  colData(cds_cm)$condition <- vcm$condition
  
  cds_cm <- cluster_cells(cds_cm, reduction_method = "UMAP")
  cds_cm <- learn_graph(cds_cm,
                        learn_graph_control = list(
                          ncenter = 65,
                          minimal_branch_len = 15,
                          prune_graph = TRUE
                        ))
  
  root_pr_node <- get_earliest_principal_node(cds_cm, "VCM3")
  
  cds_cm <- order_cells(cds_cm, root_pr_nodes = root_pr_node)
  
  p <- plot_cells(cds_cm, color_cells_by = "pseudotime",
                  label_roots = TRUE, label_leaves = TRUE,
                  label_branch_points = TRUE, graph_label_size = 5) +
    ggtitle(paste0(cond)) +
    theme(legend.text = element_text(size = 20),
          legend.title = element_text(size = 24),
          plot.title = element_text(size = 30, face = "bold", hjust = 0.5)) &
    coord_cartesian(xlim = c(-17.5, -8), ylim = c(-15, -8))
  
  cds_list[[cond]] <- cds_cm
  plot_list[[cond]] <- p
}

plot_list$W8
plot_list$W16
