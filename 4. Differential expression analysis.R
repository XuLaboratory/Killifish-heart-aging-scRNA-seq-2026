# ==============================================================================
# Differential Expression Analysis - Killifish Heart Aging Study
# 
# Purpose: To identify differentially expressed genes between W8 and W16
#          for each cell type
# 
# Input: kf_simplified.rds 
#
# Output: 
#   - DEGs_Top50.xlsx: Top 50 DEGs per cluster
#   - Volcano_*.pdf/png: Individual volcano plots
#   - Merged_DEGs.xlsx: DEGs for merged groups
#   - Volcano_Horizontal/Grid.pdf/png: Multi-panel plots
# 
# Author: Edward (Yifeng) Xu
# ==============================================================================

message("\n=== Differential Expression Analysis ===\n")

# --- Environment Setup ---

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

needed_packages <- c("Seurat", "dplyr", "ggplot2", "writexl", "readxl",
                     "ggrepel", "patchwork")

new_packages <- needed_packages[!(needed_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) BiocManager::install(new_packages)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(writexl)
  library(readxl)
  library(ggrepel)
  library(patchwork)
})

# ------------------------------------------------------------------------------
# Path Configuration
# MODIFY THESE PATHS before running:
# - Set base_dir to your data directory
# - Set out_dir to your desired output directory
# ------------------------------------------------------------------------------

# Option 1: Use current working directory (recommended)
base_dir <- getwd()

# Option 2: Specify custom path (uncomment and modify)
# base_dir <- "/path/to/your/data/"

# Output directory
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

# Define conditions
w16 <- "W16"
w8 <- "W8"

# Gene name conversion function
kf_to_human <- function(genes) {
  g <- tolower(genes)
  g <- sub("[ab]$", "", g)
  g <- toupper(g)
  return(g)
}


# --- DEG Analysis for All Clusters ---

message("\n--- Running DEG analysis for all clusters ---")

all_degs_list <- list()
target_clusters <- levels(kf)

for (cluster in target_clusters) {
  sub_obj <- subset(kf, idents = cluster)
  if (min(table(sub_obj$condition)) < 3) {
    message(paste("  Skipping", cluster, "- insufficient cells"))
    next
  }
  
  try({
    markers <- FindMarkers(sub_obj, 
                           ident.1 = w16,
                           ident.2 = w8,
                           group.by = "condition", 
                           min.pct = 0.1, 
                           logfc.threshold = 0,
                           verbose = FALSE)
    
    if(nrow(markers) > 0) {
      markers$gene <- rownames(markers)
      markers$cluster <- cluster
      markers$gene_human <- kf_to_human(markers$gene)
      all_degs_list[[cluster]] <- markers
      message(paste("  ✓", cluster, "-", nrow(markers), "DEGs"))
    }
  }, silent = TRUE)
}

# Export top 50 DEGs
if(length(all_degs_list) > 0) {
  deg_export <- lapply(all_degs_list, function(x) {
    x %>% dplyr::top_n(50, abs(avg_log2FC)) %>% 
      dplyr::select(gene, gene_human, avg_log2FC, pct.1, pct.2, p_val, p_val_adj)
  })
  write_xlsx(deg_export, paste0(out_dir, "DEGs_Top50.xlsx"))
  message("✓ Top 50 DEGs exported")
}

# --- Volcano Plot Function ---

plot_volcano <- function(deg_df, title = "Volcano Plot", 
                         fc_threshold = 0.5, pval_threshold = 0.05,
                         top_n_labels = 10, show_legend = TRUE) {
  
  plot_data <- deg_df %>%
    dplyr::mutate(
      neg_log_pval = -log10(p_val_adj + 1e-300),
      significance = dplyr::case_when(
        avg_log2FC > fc_threshold & p_val_adj < pval_threshold ~ "Up in W16",
        avg_log2FC < -fc_threshold & p_val_adj < pval_threshold ~ "Down in W16",
        TRUE ~ "NS"
      )
    )
  
  top_up <- plot_data %>%
    dplyr::filter(significance == "Up in W16") %>%
    dplyr::arrange(desc(avg_log2FC)) %>%
    head(top_n_labels)
  
  top_down <- plot_data %>%
    dplyr::filter(significance == "Down in W16") %>%
    dplyr::arrange(avg_log2FC) %>%
    head(top_n_labels)
  
  label_genes <- rbind(top_up, top_down)
  
  n_up <- sum(plot_data$significance == "Up in W16")
  n_down <- sum(plot_data$significance == "Down in W16")
  
  p <- ggplot(plot_data, aes(x = avg_log2FC, y = neg_log_pval)) +
    geom_point(aes(color = significance), alpha = 0.6, size = 1.5) +
    geom_point(data = label_genes, color = "black", size = 2) +
    geom_text_repel(
      data = label_genes,
      aes(label = gene),
      size = 3,
      max.overlaps = 20,
      box.padding = 0.5
    ) +
    geom_vline(xintercept = c(-fc_threshold, fc_threshold), 
               linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = -log10(pval_threshold), 
               linetype = "dashed", color = "gray50") +
    scale_color_manual(
      values = c("Up in W16" = "#E41A1C", 
                 "Down in W16" = "#377EB8", 
                 "NS" = "gray70"),
      name = "Regulation"
    ) +
    labs(
      title = title,
      subtitle = paste("Up:", n_up, "| Down:", n_down),
      x = "log2(Fold Change) [W16 vs W8]",
      y = "-log10(Adjusted P-value)"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 9)
    )
  
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }
  
  return(p)
}

# --- Generate Individual Volcano Plots ---

message("\n--- Generating individual volcano plots ---")

for (ct in names(all_degs_list)) {
  p <- plot_volcano(all_degs_list[[ct]], title = paste0(ct, " (W16 vs W8)"))
  ggsave(paste0(out_dir, "Volcano_", gsub("-", "_", ct), ".pdf"), 
         p, width = 8, height = 6)
  ggsave(paste0(out_dir, "Volcano_", gsub("-", "_", ct), ".png"), 
         p, width = 8, height = 6, dpi = 300)
}

message(paste("✓", length(all_degs_list), "volcano plots saved"))

# --- Merged Group Analysis (for publication) ---

message("\n--- Analyzing merged groups ---")

merged_groups <- list(
  "M"   = c("M"),
  "B"   = c("B"),
  "EC"  = c("EC"),
  "VCM" = c("VCM"),
  "Epi" = c("Epi"),
  "T"   = c("T")
)

existing_clusters <- levels(kf)
if ("EPDC" %in% existing_clusters && "Epi" %in% existing_clusters) {
  merged_groups[["Epi"]] <- c("Epi", "EPDC")
  message("  Merging Epi and EPDC")
}

merged_deg_results <- list()

for (group_name in names(merged_groups)) {
  valid_members <- intersect(merged_groups[[group_name]], existing_clusters)
  
  if (length(valid_members) == 0) next
  
  sub_obj <- subset(kf, idents = valid_members)
  
  if (min(table(sub_obj$condition)) < 5) {
    message(paste("  Skipping", group_name, "- insufficient cells"))
    next
  }
  
  tryCatch({
    markers <- FindMarkers(sub_obj, 
                           ident.1 = w16,
                           ident.2 = w8,
                           group.by = "condition", 
                           min.pct = 0.1, 
                           logfc.threshold = 0,
                           verbose = FALSE)
    
    if (nrow(markers) > 0) {
      markers$gene <- rownames(markers)
      markers$group <- group_name
      merged_deg_results[[group_name]] <- markers
      message(paste("  ✓", group_name, "-", nrow(markers), "DEGs"))
    }
  }, error = function(e) {
    message(paste("  Error:", group_name, "-", e$message))
  })
}

# Export merged results
if (length(merged_deg_results) > 0) {
  merged_export <- lapply(merged_deg_results, function(x) {
    x %>% dplyr::select(gene, avg_log2FC, pct.1, pct.2, p_val, p_val_adj) %>%
      dplyr::arrange(p_val_adj)
  })
  write_xlsx(merged_export, paste0(out_dir, "Merged_DEGs.xlsx"))
  message("✓ Merged group DEGs exported")
}

# Generate standard volcano plots for merged groups
for (group_name in names(merged_deg_results)) {
  p <- plot_volcano(merged_deg_results[[group_name]], 
                    title = paste0(group_name, " (W16 vs W8)"))
  ggsave(paste0(out_dir, "Volcano_Merged_", group_name, ".pdf"), 
         p, width = 8, height = 6)
  ggsave(paste0(out_dir, "Volcano_Merged_", group_name, ".png"), 
         p, width = 8, height = 6, dpi = 300)
}

# --- Publication-Ready Volcano Plots ---

message("\n--- Generating publication-ready plots ---")

plot_volcano_paper <- function(deg_df, title_name, 
                               fc_threshold = 0.5, pval_threshold = 0.05,
                               top_n_labels = 5,
                               show_x_label = TRUE,
                               show_y_label = TRUE,
                               is_T_cluster = FALSE) {
  
  plot_data <- deg_df %>%
    dplyr::mutate(
      neg_log_pval = -log10(p_val_adj + 1e-300),
      significance = dplyr::case_when(
        avg_log2FC > fc_threshold & p_val_adj < pval_threshold ~ "Up",
        avg_log2FC < -fc_threshold & p_val_adj < pval_threshold ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  top_up <- plot_data %>% filter(significance == "Up") %>% 
    arrange(desc(avg_log2FC)) %>% head(top_n_labels)
  top_down <- plot_data %>% filter(significance == "Down") %>% 
    arrange(avg_log2FC) %>% head(top_n_labels)
  label_genes <- rbind(top_up, top_down)
  
  n_up <- sum(plot_data$significance == "Up")
  n_down <- sum(plot_data$significance == "Down")
  
  x_label <- if(show_x_label) expression(paste("avg_log"[2], "FC")) else NULL
  y_label <- if(show_y_label) expression(paste("-log"[10], "(p_val_adj)")) else NULL
  
  if (is_T_cluster) {
    p <- ggplot(plot_data, aes(x = avg_log2FC, y = neg_log_pval)) +
      geom_point(aes(color = significance), alpha = 0.6, size = 1.2) +
      geom_text_repel(
        data = label_genes,
        aes(label = gene),
        size = 3.0,
        max.overlaps = 20,
        box.padding = 1.2,
        point.padding = 0.8,
        segment.size = 0.25,
        fontface = "italic",
        force = 7,
        direction = "y",
        min.segment.length = 0,
        max.time = 3,
        max.iter = 15000
      )
  } else {
    p <- ggplot(plot_data, aes(x = avg_log2FC, y = neg_log_pval)) +
      geom_point(aes(color = significance), alpha = 0.6, size = 1.2) +
      geom_text_repel(
        data = label_genes,
        aes(label = gene),
        size = 3.5,
        max.overlaps = 15,
        box.padding = 0.4,
        segment.size = 0.3,
        fontface = "italic"
      )
  }
  
  p <- p +
    geom_vline(xintercept = c(-fc_threshold, fc_threshold), 
               linetype = "dashed", color = "gray50", linewidth = 0.4) +
    geom_hline(yintercept = -log10(pval_threshold), 
               linetype = "dashed", color = "gray50", linewidth = 0.4) +
    scale_color_manual(
      values = c("Up" = "#E41A1C", "Down" = "#377EB8", "NS" = "gray80")
    ) +
    labs(title = title_name, x = x_label, y = y_label) +
    annotate("text", x = Inf, y = Inf, label = paste0(n_up), 
             hjust = 1.2, vjust = 1.5, size = 4.5, color = "#E41A1C", fontface = "bold") +
    annotate("text", x = -Inf, y = Inf, label = paste0(n_down), 
             hjust = -0.2, vjust = 1.5, size = 4.5, color = "#377EB8", fontface = "bold") +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
    )
  
  return(p)
}

major_cell_types <- c("M", "B", "EC", "VCM", "Epi", "T")
paper_volcano_plots <- list()
plot_order <- intersect(major_cell_types, names(merged_deg_results))

for (i in seq_along(plot_order)) {
  group_name <- plot_order[i]
  
  if (group_name %in% names(merged_deg_results)) {
    show_y <- TRUE
    
    if (group_name == "T") {
      paper_volcano_plots[[group_name]] <- plot_volcano_paper(
        merged_deg_results[[group_name]],
        title_name = group_name,
        top_n_labels = 5,
        show_x_label = TRUE,
        show_y_label = show_y,
        is_T_cluster = TRUE
      )
    } else {
      paper_volcano_plots[[group_name]] <- plot_volcano_paper(
        merged_deg_results[[group_name]],
        title_name = group_name,
        top_n_labels = 5,
        show_x_label = TRUE,
        show_y_label = show_y,
        is_T_cluster = FALSE
      )
    }
  }
}

# Multi-panel layouts
if (length(paper_volcano_plots) > 0) {
  
  # Horizontal layout
  plot_list_1x <- paper_volcano_plots
  for (i in seq_along(plot_list_1x)) {
    if (i != 1) {
      plot_list_1x[[i]] <- plot_list_1x[[i]] + 
        theme(axis.title.y = element_blank())
    }
  }
  
  combined_paper <- wrap_plots(paper_volcano_plots, nrow = 1) +
    plot_annotation(theme = theme(plot.margin = margin(5, 5, 5, 5)))
  
  ggsave(paste0(out_dir, "Volcano_Horizontal.pdf"), 
         combined_paper, width = 18, height = 3.5)
  ggsave(paste0(out_dir, "Volcano_Horizontal.png"), 
         combined_paper, width = 18, height = 3.5, dpi = 400)
  
  # Grid layout
  n_plots <- length(paper_volcano_plots)
  ncol_grid <- 3
  nrow_grid <- ceiling(n_plots / ncol_grid)
  
  plot_list_grid <- paper_volcano_plots
  for (i in seq_along(plot_list_grid)) {
    col_num <- ((i - 1) %% ncol_grid) + 1
    if (col_num != 1) {
      plot_list_grid[[i]] <- plot_list_grid[[i]] + 
        theme(axis.title.y = element_blank())
    }
  }
  
  combined_paper_grid <- wrap_plots(plot_list_grid, nrow = nrow_grid, ncol = ncol_grid) +
    plot_annotation(theme = theme(plot.margin = margin(5, 5, 5, 5)))
  
  ggsave(paste0(out_dir, "Volcano_Grid.pdf"), 
         combined_paper_grid, width = 12, height = 7)
  ggsave(paste0(out_dir, "Volcano_Grid.png"), 
         combined_paper_grid, width = 12, height = 7, dpi = 400)
  
  message("✓ Publication-ready plots saved")
}

# Export all results
if (length(all_degs_list) > 0) {
  all_export <- lapply(all_degs_list, function(x) {
    x %>% dplyr::select(gene, avg_log2FC, pct.1, pct.2, p_val, p_val_adj) %>%
      dplyr::arrange(p_val_adj)
  })
  write_xlsx(all_export, paste0(out_dir, "All_DEGs.xlsx"))
}

