suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

# =========================================================
# EQUALIZE DMSO TO PRBJN CELL COUNT AND RE-CLUSTER
#
# What this script does:
# 1. Loads a Seurat RDS
# 2. Merges sample labels into DMSO / PRBJN
# 3. Randomly downsamples DMSO to match PRBJN cell count
# 4. Saves a new equalized Seurat object as RDS
# 5. Re-runs standard Seurat clustering on the equalized object
# 6. Saves UMAP and cluster/sample plots
# =========================================================

# -------------------------
# USER INPUTS
# -------------------------
input_rds <- "path/to/seu_sample_merged_DMSO_PRBJN.rds"
out_dir   <- "equalized_dmso_prbjn_seurat"

sample_col  <- "sample"
random_seed <- 123

n_variable_features <- 3000
dims_use           <- 1:30
umap_dims          <- 1:30
cluster_resolution <- 0.8

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("Loading Seurat object...\n")
seu <- readRDS(input_rds)

if (!sample_col %in% colnames(seu@meta.data)) {
  stop(paste0("Sample column not found in metadata: ", sample_col))
}

merge_sample_labels <- function(x) {
  ifelse(
    grepl("DMSO", x, ignore.case = TRUE),
    "DMSO",
    ifelse(
      grepl("PRBJN", x, ignore.case = TRUE),
      "PRBJN",
      x
    )
  )
}

seu$sample_merged <- merge_sample_labels(as.character(seu@meta.data[[sample_col]]))

cat("\nOriginal sample labels:\n")
print(table(seu@meta.data[[sample_col]], useNA = "ifany"))

cat("\nMerged sample labels:\n")
print(table(seu$sample_merged, useNA = "ifany"))

if (!all(c("DMSO", "PRBJN") %in% unique(seu$sample_merged))) {
  stop("Both DMSO and PRBJN must be present after merging labels.")
}

set.seed(random_seed)

dmso_cells  <- colnames(seu)[seu$sample_merged == "DMSO"]
prbjn_cells <- colnames(seu)[seu$sample_merged == "PRBJN"]

n_dmso  <- length(dmso_cells)
n_prbjn <- length(prbjn_cells)

cat("\nBefore equalization:\n")
cat("DMSO cells :", n_dmso, "\n")
cat("PRBJN cells:", n_prbjn, "\n")

if (n_dmso < n_prbjn) {
  stop("DMSO has fewer cells than PRBJN. This script expects DMSO >= PRBJN.")
}

dmso_keep  <- sample(dmso_cells, size = n_prbjn, replace = FALSE)
cells_keep <- c(dmso_keep, prbjn_cells)

seu_equal <- subset(seu, cells = cells_keep)

cat("\nAfter equalization:\n")
print(table(seu_equal$sample_merged, useNA = "ifany"))

cat("\nRunning Seurat pipeline on equalized object...\n")

DefaultAssay(seu_equal) <- "RNA"

existing_layers <- tryCatch(Layers(seu_equal[["RNA"]]), error = function(e) character(0))
cat("RNA layers before processing:", paste(existing_layers, collapse = ", "), "\n")

if (!"data" %in% existing_layers) {
  seu_equal <- NormalizeData(seu_equal, verbose = FALSE)
}

seu_equal <- FindVariableFeatures(
  seu_equal,
  selection.method = "vst",
  nfeatures = n_variable_features,
  verbose = FALSE
)

seu_equal <- ScaleData(seu_equal, verbose = FALSE)
seu_equal <- RunPCA(seu_equal, npcs = max(dims_use), verbose = FALSE)
seu_equal <- FindNeighbors(seu_equal, dims = dims_use, verbose = FALSE)
seu_equal <- FindClusters(seu_equal, resolution = cluster_resolution, verbose = FALSE)
seu_equal <- RunUMAP(seu_equal, dims = umap_dims, verbose = FALSE)

seu_equal$equalized_clusters <- as.character(Idents(seu_equal))

out_rds <- file.path(out_dir, "equalized_dmso_prbjn_seurat.rds")
saveRDS(seu_equal, out_rds)

meta_out <- seu_equal@meta.data
write.csv(meta_out, file.path(out_dir, "equalized_metadata.csv"))

cluster_sample_table <- table(seu_equal$sample_merged, seu_equal$equalized_clusters)
write.csv(as.data.frame.matrix(cluster_sample_table),
          file.path(out_dir, "sample_by_cluster_table.csv"))

cluster_percent_df <- as.data.frame(cluster_sample_table)
colnames(cluster_percent_df) <- c("sample_merged", "cluster", "cell_count")

cluster_percent_df <- do.call(
  rbind,
  lapply(split(cluster_percent_df, cluster_percent_df$sample_merged), function(df) {
    df$percent_within_sample <- round(100 * df$cell_count / sum(df$cell_count), 2)
    df
  })
)

write.csv(cluster_percent_df, file.path(out_dir, "sample_by_cluster_percent.csv"), row.names = FALSE)

p1 <- DimPlot(seu_equal, group.by = "sample_merged") +
  ggtitle("Equalized object - colored by sample")

p2 <- DimPlot(seu_equal, group.by = "equalized_clusters", label = TRUE) +
  ggtitle("Equalized object - colored by clusters")

p3 <- DimPlot(seu_equal, group.by = "sample_merged", split.by = "sample_merged") +
  ggtitle("Equalized object - split by sample")

ggsave(file.path(out_dir, "umap_by_sample.png"),       p1, width = 8,  height = 6, dpi = 300)
ggsave(file.path(out_dir, "umap_by_clusters.png"),     p2, width = 8,  height = 6, dpi = 300)
ggsave(file.path(out_dir, "umap_split_by_sample.png"), p3, width = 12, height = 6, dpi = 300)

ggsave(file.path(out_dir, "umap_by_sample.pdf"),       p1, width = 8,  height = 6)
ggsave(file.path(out_dir, "umap_by_clusters.pdf"),     p2, width = 8,  height = 6)
ggsave(file.path(out_dir, "umap_split_by_sample.pdf"), p3, width = 12, height = 6)

summary_file <- file.path(out_dir, "summary.txt")
cat_lines <- c(
  "Equalized DMSO / PRBJN Seurat object summary",
  "===========================================",
  "",
  paste("Input RDS:", input_rds),
  paste("Output RDS:", out_rds),
  "",
  paste("Original DMSO cell count:", n_dmso),
  paste("Original PRBJN cell count:", n_prbjn),
  paste("Final DMSO cell count:", sum(seu_equal$sample_merged == 'DMSO')),
  paste("Final PRBJN cell count:", sum(seu_equal$sample_merged == 'PRBJN')),
  "",
  paste("Cluster resolution:", cluster_resolution),
  paste("Number of clusters:", length(unique(seu_equal$equalized_clusters)))
)
writeLines(cat_lines, summary_file)

cat("\nDone.\n")
cat("Saved files in:\n", normalizePath(out_dir), "\n")
cat("Main RDS file:\n", normalizePath(out_rds), "\n")
