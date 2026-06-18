suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

# =========================================================
# VIOLIN PLOT - ZGA Gene Expression Across Clusters
# =========================================================

# -------------------------
# USER INPUTS
# -------------------------
rds_file    <- "path/to/prbjn.rds"
out_dir     <- "ViolinPlot_ZGA_outputs"
cluster_col <- "RNA_snn_res.0.8"
assay_use   <- "RNA"

# Genes to plot
genes_of_interest <- c(
  "TPRX1", "ZSCAN4", "DUXA", "DUXB",
  "KLF17", "SLC34A2", "CCNA1", "RFPL2", "DPPA3"
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =========================================================
# LOAD SEURAT OBJECT
# =========================================================
cat("Loading Seurat object...\n")
seu <- readRDS(rds_file)
DefaultAssay(seu) <- assay_use

if (!cluster_col %in% colnames(seu@meta.data)) {
  stop(paste0("Cluster column not found in metadata: ", cluster_col))
}

Idents(seu) <- seu@meta.data[[cluster_col]]

cat("Cells :", ncol(seu), "\n")
cat("Genes :", nrow(seu), "\n")
cat("Clusters:", paste(sort(unique(as.character(Idents(seu)))), collapse = ", "), "\n")

# =========================================================
# MATCH GENES (case-insensitive)
# =========================================================
all_features <- rownames(seu)
feature_map  <- setNames(all_features, toupper(all_features))

matched <- feature_map[toupper(genes_of_interest)]
matched <- unname(matched[!is.na(matched)])
matched <- unique(matched)

missing <- genes_of_interest[!(toupper(genes_of_interest) %in% names(feature_map))]

cat("\nMatched genes :", paste(matched, collapse = ", "), "\n")
if (length(missing) > 0) {
  cat("Missing genes :", paste(missing, collapse = ", "), "\n")
}
if (length(matched) == 0) stop("No genes found in the Seurat object.")

# =========================================================
# CLUSTER COLORS (same palette as the heatmap script)
# =========================================================
cluster_levels <- sort(unique(as.character(seu@meta.data[[cluster_col]])), decreasing = FALSE)

# Numeric sort if all cluster names are numbers
suppressWarnings(num_test <- as.numeric(cluster_levels))
if (all(!is.na(num_test))) {
  cluster_levels <- cluster_levels[order(as.numeric(cluster_levels))]
}

base_colors <- c(
  "#d27d73", "#b89d22", "#68a828", "#5aa081", "#5a95d6",
  "#8579b9", "#b267a5", "#d98080", "#7fbf7b", "#80b1d3"
)

n_clusters <- length(cluster_levels)
if (n_clusters > length(base_colors)) {
  extra <- colorRampPalette(base_colors)(n_clusters)
  cluster_color_map <- setNames(extra, cluster_levels)
} else {
  cluster_color_map <- setNames(base_colors[seq_len(n_clusters)], cluster_levels)
}

# =========================================================
# THEME - clean, publication-ready
# =========================================================
violin_theme <- theme_classic(base_size = 14) +
  theme(
    axis.title.y    = element_text(face = "bold", size = 13, angle = 90),
    axis.title.x    = element_blank(),
    axis.text.x     = element_text(size = 11, angle = 45, hjust = 1, vjust = 1),
    axis.text.y     = element_text(size = 10),
    legend.position = "none",
    plot.margin     = margin(t = 2, r = 6, b = 2, l = 6, unit = "pt"),
    panel.border    = element_rect(color = "grey60", fill = NA, linewidth = 0.4),
    plot.title      = element_blank()
  )

# =========================================================
# BUILD INDIVIDUAL VIOLIN PANELS
# =========================================================
plot_list <- lapply(seq_along(matched), function(i) {
  gene    <- matched[i]
  is_last <- (i == length(matched))

  p <- VlnPlot(
    seu, features = gene, group.by = cluster_col,
    cols = cluster_color_map, pt.size = 0, sort = FALSE
  ) +
    ylab(gene) +
    violin_theme

  # Hide x-axis tick labels on all panels except the bottom one
  if (!is_last) {
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }

  p
})

# =========================================================
# STACK PANELS WITH PATCHWORK
# =========================================================
n_genes <- length(plot_list)

combined <- patchwork::wrap_plots(plot_list, ncol = 1, heights = rep(1, n_genes)) +
  patchwork::plot_annotation(
    title = "ZGA Gene Expression Across Clusters",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(b = 8))
    )
  )

# =========================================================
# SAVE
# =========================================================
panel_height <- 1.8                            # inches per gene panel
total_height <- panel_height * n_genes + 1.2   # +1.2 for title + x labels
plot_width   <- 7

out_pdf <- file.path(out_dir, "ViolinPlot_ZGA_genes_by_cluster.pdf")
out_png <- file.path(out_dir, "ViolinPlot_ZGA_genes_by_cluster.png")

ggsave(out_pdf, combined, width = plot_width, height = total_height)
ggsave(out_png, combined, width = plot_width, height = total_height, dpi = 300)

cat("\nSaved:\n", out_pdf, "\n", out_png, "\n")
cat("Done.\n")
