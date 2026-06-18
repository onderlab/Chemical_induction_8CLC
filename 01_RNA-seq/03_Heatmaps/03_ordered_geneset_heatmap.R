# Heatmap in Gene Set Order (Row Order Preserved)
# -------------------------------------------------
# 1. Loads a TPM matrix and a gene set
# 2. Keeps the gene order from the gene set (no row clustering)
# 3. Averages samples per condition
# 4. Plots a ComplexHeatmap with row z-scores

library(readxl)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(openxlsx)

# ============================================================
# 1. Read TPM matrix
# ============================================================
tpm_df_raw <- read_excel("path/to/TPM_matrix_(proteincoding).xlsx")
tpm_df     <- as.data.frame(tpm_df_raw)
rownames(tpm_df) <- tpm_df$gene_symbol
tpm_df <- tpm_df[, !colnames(tpm_df) %in% "gene_symbol"]

# ============================================================
# 2. Read gene set
# ============================================================
geneset   <- read_excel("path/to/dux4-target_geneset.xlsx")
gene_list <- as.character(geneset[[1]])

# ============================================================
# 3. Standardize gene names (uppercase, trim whitespace)
# ============================================================
rownames(tpm_df) <- toupper(trimws(rownames(tpm_df)))
gene_list        <- toupper(trimws(gene_list))

# ============================================================
# 4. Subset to gene set in order
# ============================================================
genes_in_order <- gene_list[gene_list %in% rownames(tpm_df)]
expr_subset    <- tpm_df[genes_in_order, , drop = FALSE]

# ============================================================
# 5. Define groups and assign samples
# ============================================================
groups <- c("Primed", "DMSO", "PY60_24h", "PY60_48h", "PY60_Sorted", "DUX4_OE")

sample_groups <- sapply(colnames(expr_subset), function(x) {
  grp <- NA
  for (g in groups) {
    if (grepl(g, x, ignore.case = TRUE)) {
      grp <- g
      break
    }
  }
  grp
})

# ============================================================
# 6. Compute per-group averages
# ============================================================
expr_avg <- sapply(groups, function(g) {
  idx <- which(sample_groups == g)
  if (length(idx) > 1) {
    rowMeans(expr_subset[, idx, drop = FALSE], na.rm = TRUE)
  } else if (length(idx) == 1) {
    expr_subset[, idx]
  } else {
    rep(NA, nrow(expr_subset))
  }
})

expr_avg <- as.data.frame(expr_avg)
rownames(expr_avg) <- rownames(expr_subset)
colnames(expr_avg) <- groups

# ============================================================
# 7. log2(TPM+1) -> row z-score
# ============================================================
expr_log <- t(scale(t(log2(expr_avg))))
expr_log <- expr_log[rowSums(is.na(expr_log)) == 0, ]
expr_log <- expr_log[, colSums(is.na(expr_log)) == 0]

# ============================================================
# 8. Dynamic color scale and heatmap dimensions
# ============================================================
range_vals <- range(expr_log, na.rm = TRUE)
col_fun <- colorRamp2(
  c(range_vals[1], 0, range_vals[2]),
  c("#2166AC", "white", "#B2182B")
)

heatmap_height <- unit(max(3 * nrow(expr_log), 20), "mm")  # minimum height

# ============================================================
# 9. ComplexHeatmap
# ============================================================
ht <- Heatmap(
  expr_log,
  name                = "row-scaled\nlog2(TPM+1)",
  col                 = col_fun,
  show_row_names      = FALSE,
  show_column_names   = TRUE,
  cluster_columns     = FALSE,
  cluster_rows        = FALSE,
  column_names_side   = "bottom",
  column_names_gp     = gpar(fontsize = 12),
  heatmap_width       = unit(6, "cm"),
  heatmap_height      = heatmap_height
)

# ============================================================
# 10. Export PNG and PDF
# ============================================================
png("geneset_heatmap_complex.png", width = 2000, height = 2000, res = 150)
draw(ht)
dev.off()

pdf("geneset_heatmap_complex_log2_filt.pdf", width = 18, height = 8)
draw(ht)
dev.off()

# ============================================================
# 11. Export z-score matrix
# ============================================================
write.xlsx(
  as.data.frame(expr_log),
  "geneset_zscore_matrix.xlsx",
  rowNames = TRUE
)
