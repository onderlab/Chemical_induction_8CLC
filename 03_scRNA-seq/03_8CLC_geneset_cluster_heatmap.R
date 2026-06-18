suppressPackageStartupMessages({
  library(Seurat)
  library(openxlsx)
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# =========================================================
# 8CLC GENE-SET HEATMAP ACROSS CLUSTERS
# Row-scaled average expression per cluster, with a curated
# subset of genes shown as row labels.
# =========================================================

# =========================================================
# INPUTS
# =========================================================
rds_file   <- "path/to/prbjn.rds"
excel_file <- "path/to/8CLC_genelist.xlsx"

out_pdf <- "8CLC_cluster_heatmap_res0.8_selectedlabels.pdf"
out_png <- "8CLC_cluster_heatmap_res0.8_selectedlabels.png"

cluster_col <- "RNA_snn_res.0.8"
assay_use   <- "RNA"

# =========================================================
# LOAD SEURAT
# =========================================================
cat("Loading Seurat object...\n")
seu <- readRDS(rds_file)

DefaultAssay(seu) <- assay_use

if (!cluster_col %in% colnames(seu@meta.data)) {
  stop(paste0("Cluster column not found: ", cluster_col))
}

Idents(seu) <- seu@meta.data[[cluster_col]]

cat("Dataset:", ncol(seu), "cells,", nrow(seu), "genes\n")
cat("Clusters:", paste(sort(unique(as.character(Idents(seu)))), collapse = " "), "\n")

# =========================================================
# READ GENE LIST
# =========================================================
gene_df <- openxlsx::read.xlsx(excel_file, colNames = TRUE)

gene_list <- gene_df[[1]]
gene_list <- as.character(gene_list)
gene_list <- trimws(gene_list)
gene_list <- gene_list[gene_list != "" & !is.na(gene_list)]
gene_list <- unique(gene_list)

# H3.Y -> H3Y1 symbol fix
gene_list[gene_list == "H3.Y"] <- "H3Y1"

# =========================================================
# MATCH GENES TO OBJECT
# =========================================================
all_features <- rownames(seu)
feature_map  <- setNames(all_features, toupper(all_features))

matched_genes <- feature_map[toupper(gene_list)]
matched_genes <- unname(matched_genes[!is.na(matched_genes)])
matched_genes <- unique(matched_genes)

missing_genes <- gene_list[!(toupper(gene_list) %in% names(feature_map))]

cat("Matched:", length(matched_genes), "\n")
cat("Missing:", length(missing_genes), "\n")
if (length(missing_genes) > 0) {
  print(missing_genes)
}

if (length(matched_genes) == 0) {
  stop("No matched genes found.")
}

# =========================================================
# AVERAGE EXPRESSION
# =========================================================
avg_expr <- tryCatch({
  AverageExpression(
    seu,
    assays = assay_use,
    features = matched_genes,
    group.by = cluster_col,
    slot = "data",
    verbose = FALSE
  )[[assay_use]]
}, error = function(e) {
  expr_mat <- tryCatch({
    GetAssayData(seu, assay = assay_use, slot = "data")
  }, error = function(e2) {
    GetAssayData(seu, assay = assay_use, layer = "data")
  })

  cluster_ids    <- as.character(seu@meta.data[[cluster_col]])
  cluster_levels <- sort(unique(cluster_ids))

  avg_list <- lapply(cluster_levels, function(cl) {
    cells_use <- colnames(seu)[cluster_ids == cl]
    Matrix::rowMeans(expr_mat[matched_genes, cells_use, drop = FALSE])
  })

  avg_mat <- do.call(cbind, avg_list)
  rownames(avg_mat) <- matched_genes
  colnames(avg_mat) <- cluster_levels
  avg_mat
})

# Cluster order (numeric if possible)
ordered_clusters <- colnames(avg_expr)
suppressWarnings(num_test <- as.numeric(ordered_clusters))
if (all(!is.na(num_test))) {
  ordered_clusters <- ordered_clusters[order(as.numeric(ordered_clusters))]
} else {
  ordered_clusters <- sort(ordered_clusters)
}
avg_expr <- avg_expr[, ordered_clusters, drop = FALSE]

# =========================================================
# SCALE BY GENE (row z-score, clamped)
# =========================================================
scaled_mat <- t(scale(t(as.matrix(avg_expr))))
scaled_mat[is.na(scaled_mat)] <- 0
scaled_mat[scaled_mat > 2]  <- 2
scaled_mat[scaled_mat < -1] <- -1

# =========================================================
# GENES TO LABEL (others left blank)
# =========================================================
genes_to_label <- c(
  "MBD3L2", "H3Y1", "DPPA3", "MBD3L2B", "CCNA1", "SLC34A2",
  "TRIM64B", "ALPG", "ODC1", "KHDC1L", "KLF17", "TRIM49B",
  "TRIM49C", "ZSCAN4", "MBD3L5", "SUSD2", "TRIM48", "TRIM51",
  "PRAMEF8", "TRIM49",
  "PRAMEF5", "TPRX1", "PRAMEF6", "PRAMEF12", "PRAMEF14", "TRIM60", "LEUTX",
  "PRAMEF1", "PITX1", "DUXA", "ARGFX"
)

# Only keep labels that are present in the matrix
genes_to_label <- genes_to_label[genes_to_label %in% rownames(scaled_mat)]

# Row labels: blank for everything else
row_labels <- ifelse(
  rownames(scaled_mat) %in% genes_to_label,
  rownames(scaled_mat),
  ""
)

# =========================================================
# CLUSTER COLORS
# =========================================================
cluster_levels <- colnames(scaled_mat)

cluster_colors <- c(
  "#d27d73", "#b89d22", "#68a828", "#5aa081", "#5a95d6",
  "#8579b9", "#b267a5", "#d98080", "#7fbf7b", "#80b1d3"
)

cluster_color_map <- setNames(cluster_colors[seq_along(cluster_levels)], cluster_levels)

ha <- HeatmapAnnotation(
  Clusters = factor(cluster_levels, levels = cluster_levels),
  col = list(Clusters = cluster_color_map),
  annotation_name_gp = gpar(fontsize = 18, fontface = "bold"),
  simple_anno_size = unit(0.5, "cm")
)

# =========================================================
# COLOR FUNCTION
# =========================================================
col_fun <- colorRamp2(
  c(-1, 0, 1, 2),
  c("#8e44ad", "black", "#7f7f2a", "#e6e632")
)

# =========================================================
# HEATMAP
# =========================================================
ht <- Heatmap(
  scaled_mat,
  name = "Expression",
  col = col_fun,
  top_annotation = ha,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  border = TRUE,

  row_labels = row_labels,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 16),
  row_names_max_width = unit(7, "cm"),

  column_names_gp = gpar(fontsize = 14),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 16, fontface = "bold"),
    labels_gp = gpar(fontsize = 12)
  )
)

# =========================================================
# SAVE
# =========================================================
pdf(out_pdf, width = 10, height = 16)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

png(out_png, width = 1800, height = 2600, res = 220)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

cat("Saved:\n", out_pdf, "\n", out_png, "\n")
