# ==============================================================================
# Data Annotation - Killifish Heart Aging Study
# 
# Purpose: To annotate cell clusters for downstream analyses
#          
# Input: kf_RIMD.rds
#
# Output: 
#   - kf_simplified.rds: Annotated Seurat object
#   - stacked_vlnplot.pdf: Stacked Violin plot visualizing known markers
#
# Author: Baul Yoon
# ==============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(purrr)

kf = readRDS("kf_RIMD.rds") 

# --- Gene markers search ---
kf$condition = factor(kf$condition, levels = c("W8", "W16"))
DefaultAssay(kf) <- "RNA"
kf <- JoinLayers(kf)

DimPlot(kf, split.by = "condition")

kf_markers = FindAllMarkers(kf, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
kf_top50_markers = kf_markers %>%  group_by(cluster) %>% top_n(n=50, wt = avg_log2FC) %>% arrange (cluster, desc(avg_log2FC))
library(writexl)
write_xlsx(kf_top50_markers, path = "kf_W8_W16_cluster_markers_top50.xlsx")

# --- Data Annotation ---
new.cluster.id <- c("M", "B", "EC", "remove", "B", "EC", "remove", "remove", 
                      "M", "M", "VCM", "Epi", "remove", "T", "remove", "EPDC", 
                      "Gran", "remove", "remove", "VCM", "VCM", "SMC", "Neu", "B", "remove")

names(new.cluster.id) <- levels(kf)
kf <- RenameIdents(kf, new.cluster.id)
kf <- subset(kf, idents = setdiff(unique(Idents(kf)), "remove"))

table(Idents(kf))
table(kf$condition)
saveRDS(kf, "kf_simplified.rds")

# --- Stacked Violin Plot ---
Idents(kf) <- "cluster"

modify_vlnplot <- function(obj, 
                           feature, 
                           pt.size = 0, 
                           plot.margin = unit(c(-1.2, 0, -1.2, 0), "cm"), 
                           ...) {
  VlnPlot(obj, features = feature, pt.size = pt.size, ... ) + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5), 
          axis.text.y = element_text(size = rel(0.8)), 
          plot.margin = plot.margin) 
}

extract_max <- function(p){
  ymax <- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}

StackedVlnPlot <- function(obj, features, pt.size = 0, ...) {
  
  plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj, feature = x, pt.size = pt.size, ...))
  
  plot_list[[length(plot_list)]] <- plot_list[[length(plot_list)]] +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          axis.ticks.x = element_line())
  
  ymaxs <- purrr::map_dbl(plot_list, extract_max)
  plot_list <- purrr::map2(plot_list, ymaxs, function(x, y) {
    x + scale_y_continuous(breaks = c(y)) + expand_limits(y = y)
  })
  
  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1) & 
    theme(plot.spacing = unit(0, "pt"))
  
  return(p)
}

feature <- c("tnfa", "mrc1b", "cd79a", "lyve1b", "nppb", "tbx18", "lck",
            "meis2a", "gfi1aa", "elnb", "atp1a3a")
final_plot <- StackedVlnPlot(obj = kf, features = feature)
ggsave("stacked_vlnplot.pdf", plot = final_plot, width = 10, height = 20)
