suppressPackageStartupMessages({
  library(Seurat)
  library(openxlsx)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
  library(tibble)
})

# =========================================================
# MULTI-GENE UMAP / FEATUREPLOT / DOTPLOT / MODULE SCORES
# For an 8CLC gene set:
#   - cluster UMAP
#   - FeaturePlots (core + selected genes)
#   - module scores (core / selected / full set)
#   - DotPlot of a curated marker panel
#   - violin / ridge plots
#   - average-expression tables
# =========================================================

# =========================================================
# USER INPUTS
# =========================================================
rds_file   <- "path/to/prbjn.rds"
excel_file <- "path/to/8CLC_genelist.xlsx"
out_dir    <- "8CLC_UMAP_outputs"

cluster_col          <- "RNA_snn_res.0.8"
assay_use            <- "RNA"
cluster_to_highlight <- "8"   # change if needed

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =========================================================
# HELPER FUNCTIONS
# =========================================================
square_umap_theme <- theme_classic(base_size = 18) +
  theme(
    plot.title   = element_text(face = "bold", size = 20, hjust = 0.5),
    axis.title   = element_text(size = 18),
    axis.text    = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 14)
  )

save_square_plot <- function(plot_obj, file_prefix, width = 6.5, height = 6) {
  ggsave(paste0(file_prefix, ".pdf"), plot = plot_obj, width = width, height = height)
  ggsave(paste0(file_prefix, ".png"), plot = plot_obj, width = width, height = height, dpi = 300)
}

get_expr_matrix <- function(seu, assay_use = "RNA") {
  tryCatch({
    GetAssayData(seu, assay = assay_use, slot = "data")
  }, error = function(e1) {
    tryCatch({
      GetAssayData(seu, assay = assay_use, layer = "data")
    }, error = function(e2) {
      stop("Could not access normalized expression matrix from assay/layer 'data'.")
    })
  })
}

get_average_expression_manual <- function(seu, features, cluster_col, assay_use = "RNA") {
  expr_mat       <- get_expr_matrix(seu, assay_use = assay_use)
  cluster_ids    <- as.character(seu@meta.data[[cluster_col]])
  cluster_levels <- sort(unique(cluster_ids))

  avg_list <- lapply(cluster_levels, function(cl) {
    cells_use <- colnames(seu)[cluster_ids == cl]
    if (length(cells_use) == 1) {
      as.numeric(expr_mat[features, cells_use, drop = FALSE])
    } else {
      Matrix::rowMeans(expr_mat[features, cells_use, drop = FALSE])
    }
  })

  avg_mat <- do.call(cbind, avg_list)
  rownames(avg_mat) <- features
  colnames(avg_mat) <- cluster_levels
  avg_mat
}

# =========================================================
# LOAD OBJECT
# =========================================================
cat("Loading Seurat object...\n")
seu <- readRDS(rds_file)
DefaultAssay(seu) <- assay_use

if (!cluster_col %in% colnames(seu@meta.data)) {
  stop(paste0("Cluster column not found in metadata: ", cluster_col))
}

if (!"umap" %in% names(seu@reductions)) {
  stop("No UMAP reduction found in Seurat object.")
}

Idents(seu) <- seu@meta.data[[cluster_col]]

cat("Cells:", ncol(seu), "\n")
cat("Genes:", nrow(seu), "\n")
cat("Clusters:", paste(sort(unique(as.character(Idents(seu)))), collapse = ", "), "\n")

# =========================================================
# READ EXCEL GENE LIST
# =========================================================
cat("Reading gene list from Excel...\n")
gene_df <- openxlsx::read.xlsx(excel_file, colNames = TRUE)

gene_list <- gene_df[[1]]
gene_list <- as.character(gene_list)
gene_list <- trimws(gene_list)
gene_list <- gene_list[gene_list != "" & !is.na(gene_list)]
gene_list <- unique(gene_list)

# common symbol fix
gene_list[gene_list == "H3.Y"] <- "H3Y1"

all_features <- rownames(seu)
feature_map  <- setNames(all_features, toupper(all_features))

matched_genes <- feature_map[toupper(gene_list)]
matched_genes <- unname(matched_genes[!is.na(matched_genes)])
matched_genes <- unique(matched_genes)

missing_genes <- gene_list[!(toupper(gene_list) %in% names(feature_map))]

cat("Matched genes:", length(matched_genes), "\n")
cat("Missing genes:", length(missing_genes), "\n")
if (length(missing_genes) > 0) {
  cat("Missing genes:\n")
  print(missing_genes)
}

if (length(matched_genes) == 0) {
  stop("No genes from Excel file were found in the Seurat object.")
}

write.table(
  data.frame(Missing_Genes = missing_genes),
  file = file.path(out_dir, "missing_genes_from_excel.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# =========================================================
# DEFINE GENE SETS
# =========================================================
core_genes <- c("TPRX1", "LEUTX", "ZSCAN4", "DUXA", "KHDC1L", "DPPA3")
core_genes <- core_genes[core_genes %in% rownames(seu)]

selected_genes <- c(
  "MBD3L2", "H3Y1", "DPPA3", "MBD3L2B", "CCNA1", "SLC34A2", "NLRP2",
  "TRIM64B", "ALPG", "ODC1", "KHDC1L", "KLF17", "TRIM49B",
  "TRIM49C", "ZSCAN4", "MBD3L5", "SUSD2", "TRIM48", "TRIM51",
  "PRAMEF8", "TRIM49",
  "PRAMEF5", "TPRX1", "PRAMEF6", "PRAMEF12", "PRAMEF14", "TRIM60", "LEUTX",
  "PRAMEF1", "PITX1", "DUXA", "ARGFX", "TFAP2A", "TFAP2B", "TFAP2C", "GATA3",
  "OTX2", "CTCFL", "FOSL2", "JUNB", "FRA2", "FRA1", "BATF", "EIF1AX", "EIF1AY",
  "SHROOM3", "RFPL4B", "PRRT2", "FASN", "EPHB4", "PTMA", "PCYT2", "MYCBP",
  "INPP5A", "TLE1", "TLE2", "TLE3", "TLE4", "RNF4", "CTCFL", "ANXA2", "ANXA1", "PBX1"
)
selected_genes <- selected_genes[selected_genes %in% rownames(seu)]

dotplot_genes_custom <- c(
  "CCNA1", "TPRX1", "KDM4E", "LEUTX", "DUXA", "DUXB", "KHDC1L",
  "MBD3L2", "PRAMEF1", "DPPA3", "KLF17",
  "RFPL2", "RFPL4A", "RFPL4B", "SLC34A2",
  "TRIM43", "TRIM43B",
  "ZNF296", "ZSCAN4"
)
dot_features <- dotplot_genes_custom[dotplot_genes_custom %in% rownames(seu)]

cat("Core genes found:", paste(core_genes, collapse = ", "), "\n")
cat("Selected genes found:", length(selected_genes), "\n")
cat("DotPlot genes found:", paste(dot_features, collapse = ", "), "\n")

# =========================================================
# MODULE SCORES
# =========================================================
if (length(core_genes) >= 2) {
  seu <- AddModuleScore(seu, features = list(core_genes), name = "Core8CLC_Score")
}
if (length(selected_genes) >= 2) {
  seu <- AddModuleScore(seu, features = list(selected_genes), name = "Selected8CLC_Score")
}
if (length(matched_genes) >= 2) {
  seu <- AddModuleScore(seu, features = list(matched_genes), name = "Full8CLC_Score")
}

# =========================================================
# DIMPLOT: CLUSTERS
# =========================================================
p_cluster <- DimPlot(
  seu, reduction = "umap", group.by = cluster_col, label = TRUE, repel = TRUE
) +
  ggtitle(paste0("UMAP - ", cluster_col)) +
  coord_fixed() +
  square_umap_theme

save_square_plot(p_cluster, file.path(out_dir, "UMAP_clusters_square"), width = 7, height = 6)

# =========================================================
# FEATUREPLOT: CORE GENES
# =========================================================
if (length(core_genes) > 0) {
  p_core <- FeaturePlot(
    seu, features = core_genes, reduction = "umap",
    cols = c("lightgrey", "red"), pt.size = 0.35,
    order = TRUE, combine = TRUE, ncol = 2
  ) &
    coord_fixed() & square_umap_theme

  ggsave(file.path(out_dir, "UMAP_core_genes_featureplot_square.pdf"), p_core, width = 11, height = 10)
  ggsave(file.path(out_dir, "UMAP_core_genes_featureplot_square.png"), p_core, width = 11, height = 10, dpi = 300)

  for (g in core_genes) {
    p_single <- FeaturePlot(
      seu, features = g, reduction = "umap",
      cols = c("lightgrey", "red"), pt.size = 0.4, order = TRUE
    ) +
      coord_fixed() + square_umap_theme

    save_square_plot(p_single, file.path(out_dir, paste0("UMAP_", g, "_square")), width = 6.5, height = 6)
  }
}

# =========================================================
# FEATUREPLOT: SELECTED GENES (in chunks of 6)
# =========================================================
if (length(selected_genes) > 0) {
  gene_chunks <- split(selected_genes, ceiling(seq_along(selected_genes) / 6))

  for (i in seq_along(gene_chunks)) {
    genes_now <- gene_chunks[[i]]

    p_sel <- FeaturePlot(
      seu, features = genes_now, reduction = "umap",
      cols = c("lightgrey", "red"), pt.size = 0.35,
      order = TRUE, ncol = 2, combine = TRUE
    ) &
      coord_fixed() & square_umap_theme

    ggsave(file.path(out_dir, paste0("UMAP_selected_genes_part", i, "_square.pdf")), p_sel, width = 11, height = 10)
    ggsave(file.path(out_dir, paste0("UMAP_selected_genes_part", i, "_square.png")), p_sel, width = 11, height = 10, dpi = 300)
  }
}

# =========================================================
# FEATUREPLOT: MODULE SCORES
# =========================================================
score_features <- c()
if ("Core8CLC_Score1"     %in% colnames(seu@meta.data)) score_features <- c(score_features, "Core8CLC_Score1")
if ("Selected8CLC_Score1" %in% colnames(seu@meta.data)) score_features <- c(score_features, "Selected8CLC_Score1")
if ("Full8CLC_Score1"     %in% colnames(seu@meta.data)) score_features <- c(score_features, "Full8CLC_Score1")

if (length(score_features) > 0) {
  p_scores <- FeaturePlot(
    seu, features = score_features, reduction = "umap",
    cols = c("lightgrey", "blue", "yellow"), pt.size = 0.35,
    order = TRUE, combine = TRUE, ncol = 2
  ) &
    coord_fixed() & square_umap_theme

  ggsave(file.path(out_dir, "UMAP_module_scores_square.pdf"), p_scores, width = 12, height = 10)
  ggsave(file.path(out_dir, "UMAP_module_scores_square.png"), p_scores, width = 12, height = 10, dpi = 300)

  for (feat in score_features) {
    p_single <- FeaturePlot(
      seu, features = feat, reduction = "umap",
      cols = c("lightgrey", "blue", "yellow"), pt.size = 0.4, order = TRUE
    ) +
      coord_fixed() + square_umap_theme

    save_square_plot(p_single, file.path(out_dir, paste0(feat, "_square")), width = 6.5, height = 6)
  }
}

# =========================================================
# DOTPLOT: CUSTOM GENE PANEL
# =========================================================
if (length(dot_features) > 0) {
  p_dot <- DotPlot(seu, features = dot_features, group.by = cluster_col) +
    RotatedAxis() +
    theme_bw(base_size = 14) +
    theme(
      plot.title  = element_text(face = "bold", hjust = 0.5, size = 18),
      axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 14),
      axis.title  = element_text(size = 16),
      legend.title = element_text(size = 14),
      legend.text  = element_text(size = 12)
    ) +
    ggtitle("DotPlot of selected 8CLC markers across clusters")

  ggsave(file.path(out_dir, "DotPlot_8CLC_markers_custom.pdf"), p_dot, width = 11, height = 6)
  ggsave(file.path(out_dir, "DotPlot_8CLC_markers_custom.png"), p_dot, width = 11, height = 6, dpi = 300)
}

# =========================================================
# VIOLIN PLOTS: MODULE SCORES
# =========================================================
if (length(score_features) > 0) {
  p_vln <- VlnPlot(seu, features = score_features, group.by = cluster_col, pt.size = 0) &
    theme_classic(base_size = 16) &
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text  = element_text(size = 12),
      axis.title = element_text(size = 14)
    )

  ggsave(file.path(out_dir, "Violin_module_scores.pdf"), p_vln, width = 4 * length(score_features), height = 5)
  ggsave(file.path(out_dir, "Violin_module_scores.png"), p_vln, width = 4 * length(score_features), height = 5, dpi = 300)
}

# =========================================================
# RIDGEPLOTS: CORE GENES
# =========================================================
if (length(core_genes) > 0) {
  ridge_features <- core_genes[1:min(length(core_genes), 6)]

  p_ridge <- RidgePlot(seu, features = ridge_features, group.by = cluster_col, ncol = 2) &
    theme_classic(base_size = 16)

  ggsave(file.path(out_dir, "RidgePlot_core_genes.pdf"), p_ridge,
         width = 10, height = 3 * ceiling(length(ridge_features) / 2))
  ggsave(file.path(out_dir, "RidgePlot_core_genes.png"), p_ridge,
         width = 10, height = 3 * ceiling(length(ridge_features) / 2), dpi = 300)
}

# =========================================================
# HIGHLIGHT CLUSTER
# =========================================================
if (cluster_to_highlight %in% unique(as.character(seu@meta.data[[cluster_col]]))) {
  cells_hi <- colnames(seu)[as.character(seu@meta.data[[cluster_col]]) == cluster_to_highlight]

  p_hi <- DimPlot(
    seu, reduction = "umap", cells.highlight = cells_hi,
    cols.highlight = "red", cols = "grey80", sizes.highlight = 0.5
  ) +
    ggtitle(paste0("Highlighted cluster ", cluster_to_highlight)) +
    coord_fixed() + square_umap_theme

  save_square_plot(p_hi, file.path(out_dir, paste0("UMAP_highlight_cluster_", cluster_to_highlight, "_square")),
                   width = 7, height = 6)
}

# =========================================================
# AVERAGE EXPRESSION TABLES
# =========================================================
avg_features_general <- unique(c(core_genes, selected_genes))
if (length(avg_features_general) > 0) {
  avg_expr_general <- tryCatch({
    AverageExpression(seu, assays = assay_use, features = avg_features_general,
                      group.by = cluster_col, slot = "data", verbose = FALSE)[[assay_use]]
  }, error = function(e) {
    get_average_expression_manual(seu, avg_features_general, cluster_col, assay_use)
  })

  avg_df_general <- tibble::rownames_to_column(as.data.frame(avg_expr_general), "Gene")
  write.csv(avg_df_general, file.path(out_dir, "AverageExpression_selected_genes_by_cluster.csv"), row.names = FALSE)
}

if (length(dot_features) > 0) {
  avg_dot <- tryCatch({
    AverageExpression(seu, assays = assay_use, features = dot_features,
                      group.by = cluster_col, slot = "data", verbose = FALSE)[[assay_use]]
  }, error = function(e) {
    get_average_expression_manual(seu, dot_features, cluster_col, assay_use)
  })

  avg_dot_df <- tibble::rownames_to_column(as.data.frame(avg_dot), "Gene")
  write.csv(avg_dot_df, file.path(out_dir, "AverageExpression_dotplot_custom_genes.csv"), row.names = FALSE)
}

# =========================================================
# SAVE UPDATED SEURAT OBJECT
# =========================================================
saveRDS(seu, file.path(out_dir, "seu_with_8CLC_scores.rds"))

cat("\nAll outputs saved in:\n", normalizePath(out_dir), "\n")
