# ==============================================================================
# Cell Type Proportion Analysis - Killifish Heart Aging Study
# 
# Purpose: To analyze cell type composition changes during cardiac aging 
#          (W8 vs W16)
# 
# Input: kf_simplified.rds
# Output: 
#   - Individual_*.png: Per-cluster comparison plots
#   - Combined_5col.png: 5-column layout
#   - Cell_Proportions.xlsx: Statistical data table
# 
# Author: Edward (Yifeng) Xu
# ==============================================================================

message("\n=== Cell Proportion Analysis ===\n")

# --- Environment Setup ---

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

needed_packages <- c("Seurat", "dplyr", "ggplot2", "writexl", 
                     "RColorBrewer", "patchwork")

new_packages <- needed_packages[!(needed_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) BiocManager::install(new_packages)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(writexl)
  library(RColorBrewer)
  library(patchwork)
})

base_dir <- "/path/to/your/data/"
out_dir <- file.path(base_dir, "results/")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# --- Data Loading ---

rds_path <- file.path(base_dir, "kf_simplified.rds")
if(file.exists(rds_path)) {
  kf <- readRDS(rds_path)
  message("✓ Data loaded successfully")
  message(paste("  Total cells:", ncol(kf)))
  message(paste("  Cell types:", paste(levels(Idents(kf)), collapse = ", ")))
} else {
  stop(paste("Error: File not found at", rds_path))
}

# --- Cell Proportion Calculation ---

cell_stats <- table(Idents(kf), kf$condition) %>% 
  as.data.frame() %>%
  dplyr::rename(Cluster = Var1, Condition = Var2, Count = Freq) %>%
  dplyr::mutate(Condition = factor(Condition, levels = c("W8", "W16"))) %>% 
  dplyr::group_by(Condition) %>%
  dplyr::mutate(Percentage = Count / sum(Count) * 100)

message("\n✓ Cell proportions calculated")

# --- Color Palette Definition ---
cluster_colors <- c(
  "VCM" = "#E41A1C",
  "EC" = "#377EB8",
  "Epi" = "#4DAF4A",
  "EPDC" = "#B2DF8A",
  "M" = "#FF7F00",
  "B" = "#A65628",
  "T" = "#66C2A5",
  "Gran" = "#8DA0CB",
  "SMC" = "#E5C494",
  "Neu" = "#B3B3B3"
)

all_clusters <- unique(cell_stats$Cluster)
colors_final <- cluster_colors[as.character(all_clusters)]
missing_colors <- is.na(colors_final)
if(any(missing_colors)) {
  n_missing <- sum(missing_colors)
  extra_colors <- colorRampPalette(brewer.pal(8, "Set3"))(n_missing)
  colors_final[missing_colors] <- extra_colors
}
names(colors_final) <- all_clusters

condition_colors <- c("W8" = "#4292C6", "W16" = "#EF3B2C")

# --- Visualization 2: Individual Comparison Plots ---

plot_list_individual <- list()

for (cluster_name in all_clusters) {
  cluster_data <- cell_stats %>% 
    dplyr::filter(Cluster == cluster_name)
  
  if (nrow(cluster_data) == 0) next
  
  max_percentage <- max(cluster_data$Percentage)
  y_max <- ceiling(max_percentage / 5) * 5 + 5
  y_breaks <- seq(0, y_max, by = 5)
  
  p_single <- ggplot(cluster_data, aes(x = Condition, y = Percentage, fill = Condition)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", 
             linewidth = 0.4, width = 0.7) +
    geom_text(aes(label = sprintf("%.2f", Percentage)), 
              vjust = -0.5, size = 5, fontface = "bold", color = "black") +
    scale_fill_manual(values = condition_colors, name = "") +
    scale_y_continuous(breaks = y_breaks, expand = c(0, 0), limits = c(0, y_max)) +
    labs(title = cluster_name, y = "Percentage [%]", x = "") +
    theme_classic(base_size = 16) +
    theme(
      plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
      axis.title.y = element_text(size = 18, face = "bold"),
      axis.text.x = element_text(size = 20, face = "bold", color = "black"),
      axis.text.y = element_text(size = 14, color = "black"),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.major.x = element_blank(),
      legend.position = "none",
      plot.margin = margin(10, 10, 10, 10)
    )
  
  ggsave(paste0(out_dir, "Individual_", cluster_name, ".png"), 
         plot = p_single, width = 4.5, height = 4.5, dpi = 300)
  
  plot_list_individual[[cluster_name]] <- p_single
}

message(paste("✓", length(plot_list_individual), "individual plots saved"))

# --- Visualization 3: Combined Layouts (sorted by proportion) ---

if (length(plot_list_individual) > 0) {
  
  cluster_total_prop <- cell_stats %>%
    dplyr::group_by(Cluster) %>%
    dplyr::summarise(total_prop = sum(Percentage)) %>%
    dplyr::arrange(desc(total_prop))
  
  cluster_order_sorted <- as.character(cluster_total_prop$Cluster)
  plot_list_sorted <- plot_list_individual[cluster_order_sorted]
  
  # Layout: 5 columns
  n_plots <- length(plot_list_sorted)
  ncol_layout1 <- 5
  nrow_layout1 <- ceiling(n_plots / ncol_layout1)
  
  plot_list_5x <- plot_list_sorted
  for (i in seq_along(plot_list_5x)) {
    col_num <- ((i - 1) %% ncol_layout1) + 1
    if (col_num != 1) {
      plot_list_5x[[i]] <- plot_list_5x[[i]] + 
        theme(axis.title.y = element_blank())
    }
  }
  
  combined_layout1 <- wrap_plots(plot_list_5x, nrow = nrow_layout1, ncol = ncol_layout1) +
    plot_annotation(
      title = "Cell Type Proportions: W8 vs W16",
      theme = theme(
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        plot.margin = margin(10, 10, 10, 10)
      )
    )
  
  ggsave(paste0(out_dir, "Combined_5col.png"), 
         combined_layout1, width = 4 * ncol_layout1, height = 4 * nrow_layout1, dpi = 300)
  
  message("✓ Combined layouts saved")
}

# --- Export Statistical Data ---

write_xlsx(cell_stats, paste0(out_dir, "Cell_Proportions.xlsx"))
message("✓ Statistical data exported")





