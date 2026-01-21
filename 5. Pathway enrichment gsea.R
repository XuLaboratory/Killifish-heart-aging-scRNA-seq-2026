# ==============================================================================
# Pathway Enrichment Analysis (GSEA) - Killifish Heart Aging Study
# 
# Function: Performs Gene Set Enrichment Analysis using human orthologs,
#           generates ridge plots showing pathway distributions
# 
# Input: kf_simplified.rds (Seurat object with merged cell types, RBC removed)
# Output: 
#   - Q3_GSEA_Results.xlsx: Complete GSEA results for all groups
#   - Q3_GSEA_RidgePlots.pdf: Ridge plots showing enrichment distributions
#   - Q3_GO_Gene_Lists.xlsx: Core enrichment genes (top 10 per pathway)
#   - Q3_RidgePlot_Gene_Summary.xlsx: Top 5 genes per pathway for ridge plots
# 
# Author: Edward (Yifeng Xu)
# ==============================================================================

message("\n=== GSEA Pathway Enrichment Analysis ===\n")

# --- Environment Setup ---

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

needed_packages <- c("Seurat", "dplyr", "ggplot2", "writexl",
                     "clusterProfiler", "org.Hs.eg.db", "enrichplot",
                     "tidyr")

new_packages <- needed_packages[!(needed_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) BiocManager::install(new_packages)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(writexl)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(tidyr)
  library(future)
})

plan("sequential")

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
cond_treatment <- "W16"
cond_control <- "W8"

# Gene name conversion function
kf_to_human <- function(genes) {
  g <- tolower(genes)
  g <- sub("[ab]$", "", g)
  g <- toupper(g)
  return(g)
}

# --- Define Broad Groups for Pathway Analysis# ---

broad_groups_pathway <- list(
  "VCM_All"    = c("VCM"),
  "EC_All"     = c("EC"),
  "Epi_All"    = c("Epi", "EPDC"),
  "Immune_All" = c("M", "T", "B")
)

# Check which groups exist
valid_groups <- list()
for (group_name in names(broad_groups_pathway)) {
  members <- broad_groups_pathway[[group_name]]
  valid_members <- intersect(members, levels(kf))
  if (length(valid_members) > 0) {
    valid_groups[[group_name]] <- valid_members
  }
}

display_names <- c(
  "VCM_All"    = "VCM",
  "EC_All"     = "EC",
  "Epi_All"    = "Epi",
  "Immune_All" = "Immune"
)

message(paste("✓ Analyzing", length(valid_groups), "cell type groups"))

# --- GSEA Analysis for Each Group ---

gsea_results_list <- list()
plot_list_ridge <- list()
gsea_objects_list <- list()
gene_mapping_list <- list()

for (group_name in names(valid_groups)) {
  
  valid_members <- valid_groups[[group_name]]
  display_name <- display_names[group_name]
  
  message(paste("\n--- Analyzing", display_name, "---"))
  
  sub_obj <- subset(kf, idents = valid_members)
  
  if (min(table(sub_obj$condition)) < 5) {
    message("  Insufficient cells, skipping")
    next
  }
  
  tryCatch({
    # Step 1: Differential expression
    markers <- FindMarkers(sub_obj, 
                           ident.1 = cond_treatment,
                           ident.2 = cond_control,
                           group.by = "condition", 
                           logfc.threshold = 0,
                           min.pct = 0.05,
                           verbose = FALSE)
    
    # Step 2: Gene name conversion
    markers$gene_kf <- rownames(markers)
    markers$gene_human <- kf_to_human(markers$gene_kf)
    markers_clean <- markers[markers$gene_human != "", ]
    
    # Step 3: Map to ENTREZID
    gene_map <- bitr(markers_clean$gene_human, 
                     fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db,
                     drop = TRUE)
    markers_clean <- merge(markers_clean, gene_map, by.x = "gene_human", by.y = "SYMBOL")
    
    # Handle duplicates
    markers_unique <- markers_clean %>% 
      group_by(ENTREZID) %>% 
      filter(abs(avg_log2FC) == max(abs(avg_log2FC))) %>% 
      ungroup()
    
    # Save gene mapping
    gene_mapping_list[[group_name]] <- markers_unique %>%
      dplyr::select(gene_kf, gene_human, ENTREZID, avg_log2FC, p_val_adj)
    
    message(paste("  ✓ Mapped", nrow(markers_unique), "genes to ENTREZID"))
    
    # Step 4: Create ranked gene list
    gene_list <- markers_unique$avg_log2FC
    names(gene_list) <- markers_unique$ENTREZID
    gene_list <- sort(gene_list, decreasing = TRUE)
    
    # Step 5: Run GSEA
    gsea_res <- gseGO(
      geneList     = gene_list,
      OrgDb        = org.Hs.eg.db,
      ont          = "BP",
      minGSSize    = 10,
      maxGSSize    = 500,
      pvalueCutoff = 0.05,
      verbose      = FALSE,
      seed         = 123
    )
    
    if (!is.null(gsea_res) && nrow(gsea_res) > 0) {
      gsea_objects_list[[group_name]] <- gsea_res
      
      # Ridge plot
      p_ridge <- ridgeplot(gsea_res, showCategory = 10) + 
        ggtitle(paste0(display_name, ": Enrichment Distribution")) +
        geom_vline(xintercept = 0, linetype="dashed", color="red") + 
        theme(
          plot.title = element_text(size = 20, face = "bold"),
          axis.text.y = element_text(size = 19, face = "bold"),
          axis.title.x = element_text(size = 14, face = "bold"),
          axis.text.x = element_text(size = 16),
          legend.title = element_text(size = 16, face = "bold"),
          legend.text = element_text(size = 14),
          legend.key.size = unit(1.0, "cm")
        )
      
      plot_list_ridge[[group_name]] <- p_ridge
      
      res_df <- as.data.frame(gsea_res)
      res_df$Group <- display_name
      gsea_results_list[[group_name]] <- res_df
      
      message(paste("  ✓ Found", nrow(gsea_res), "significant pathways"))
    } else {
      message("  No significant pathways")
    }
    
  }, error = function(e) { 
    message(paste("  Error:", e$message))
  })
}

# --- Save GSEA Results ---

if (length(gsea_results_list) > 0) {
  write_xlsx(do.call(rbind, gsea_results_list), paste0(out_dir, "Q3_GSEA_Results.xlsx"))
  message("\n✓ GSEA results exported")
  
  pdf(paste0(out_dir, "Q3_GSEA_RidgePlots.pdf"), width = 10, height = 8)
  for (n in names(plot_list_ridge)) print(plot_list_ridge[[n]])
  dev.off()
  message("✓ Ridge plots saved")
}

# --- Extract Core Enrichment Genes ---

message("\n--- Extracting core enrichment genes ---")

extract_go_genes <- function(gsea_obj, gene_mapping_df, group_name, top_n = 10) {
  
  if (is.null(gsea_obj) || nrow(gsea_obj) == 0) return(NULL)
  
  gsea_df <- as.data.frame(gsea_obj)
  go_gene_list <- list()
  
  for (i in 1:nrow(gsea_df)) {
    go_id <- gsea_df$ID[i]
    go_description <- gsea_df$Description[i]
    core_genes <- gsea_df$core_enrichment[i]
    nes <- gsea_df$NES[i]
    pvalue <- gsea_df$pvalue[i]
    p_adjust <- gsea_df$p.adjust[i]
    
    gene_ids <- strsplit(core_genes, "/")[[1]]
    
    gene_info <- gene_mapping_df %>%
      dplyr::filter(ENTREZID %in% gene_ids) %>%
      dplyr::arrange(desc(abs(avg_log2FC))) %>%
      dplyr::mutate(
        GO_ID = go_id,
        GO_Description = go_description,
        Group = group_name,
        NES = nes,
        pvalue = pvalue,
        p_adjust = p_adjust
      ) %>%
      dplyr::select(GO_ID, GO_Description, Group, NES, pvalue, p_adjust,
                    gene_kf, gene_human, ENTREZID, avg_log2FC, p_val_adj)
    
    if (nrow(gene_info) > top_n) {
      gene_info <- gene_info %>% head(top_n)
    }
    
    go_gene_list[[go_id]] <- gene_info
  }
  
  all_genes_df <- do.call(rbind, go_gene_list)
  return(all_genes_df)
}

all_go_genes <- list()

for (group_name in names(gsea_objects_list)) {
  
  display_name <- display_names[group_name]
  
  go_genes <- extract_go_genes(
    gsea_obj = gsea_objects_list[[group_name]],
    gene_mapping_df = gene_mapping_list[[group_name]],
    group_name = display_name,
    top_n = 10
  )
  
  if (!is.null(go_genes)) {
    all_go_genes[[display_name]] <- go_genes
    message(paste("  ✓", display_name, "-", nrow(go_genes), "gene entries"))
  }
}

if (length(all_go_genes) > 0) {
  write_xlsx(all_go_genes, paste0(out_dir, "Q3_GO_Gene_Lists.xlsx"))
  message("✓ Core enrichment genes exported")
}

# --- Ridge Plot Gene Summary (Top 5 per pathway) ---

message("\n--- Creating ridge plot gene summaries ---")

ridge_gene_summary <- list()

for (group_name in names(gsea_objects_list)) {
  
  display_name <- display_names[group_name]
  
  gsea_obj <- gsea_objects_list[[group_name]]
  gene_map <- gene_mapping_list[[group_name]]
  
  if (is.null(gsea_obj) || nrow(gsea_obj) == 0) next
  
  gsea_df <- as.data.frame(gsea_obj) %>% head(10)
  summary_list <- list()
  
  for (i in 1:nrow(gsea_df)) {
    go_desc <- gsea_df$Description[i]
    core_genes <- gsea_df$core_enrichment[i]
    gene_ids <- strsplit(core_genes, "/")[[1]]
    
    gene_details <- gene_map %>%
      dplyr::filter(ENTREZID %in% gene_ids) %>%
      dplyr::arrange(desc(abs(avg_log2FC))) %>%
      head(5) %>%
      dplyr::mutate(
        Pathway = go_desc,
        Rank = row_number()
      ) %>%
      dplyr::select(Pathway, Rank, gene_kf, gene_human, avg_log2FC, p_val_adj)
    
    summary_list[[i]] <- gene_details
  }
  
  ridge_gene_summary[[display_name]] <- do.call(rbind, summary_list)
}

if (length(ridge_gene_summary) > 0) {
  write_xlsx(ridge_gene_summary, 
             paste0(out_dir, "Q3_RidgePlot_Gene_Summary.xlsx"))
  message("✓ Ridge plot gene summaries exported")
}

# ---# Summary ---

message("\n=== Analysis Complete ===")
message(paste("Output directory:", out_dir))
message("Generated files:")
message("  - Q3_GSEA_Results.xlsx")
message("  - Q3_GSEA_RidgePlots.pdf")
message("  - Q3_GO_Gene_Lists.xlsx")
message("  - Q3_RidgePlot_Gene_Summary.xlsx")
message(paste("\nAnalyzed", length(gsea_results_list), "cell type groups"))
