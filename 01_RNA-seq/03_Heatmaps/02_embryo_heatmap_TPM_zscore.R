# Embryo Stage Heatmap from TPM (Group-Averaged, Z-score Scaled)
# ---------------------------------------------------------------
# 1. Reads a TPM matrix merged with the Yan et al. 2013 dataset
# 2. Subsets to a user-defined gene set (in dataset order)
# 3. Averages samples per developmental group
# 4. Plots a clustered/non-clustered heatmap of row z-scores

library(readxl)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(pheatmap)

# ============================================================
# 1. Read the TPM matrix
# ============================================================
tpm_full <- read_excel(
  "path/to/Merged_yan_et_al_2013_dataset_and_this_study_TPM_matrix.xlsx"
)

gene_symbols <- tpm_full[[1]]
tpm_df <- as.data.frame(tpm_full[, -1])
rownames(tpm_df) <- toupper(trimws(as.character(gene_symbols)))

# ============================================================
# 2. Load gene set
# ============================================================
geneset   <- read_excel("path/to/zygote-2cell_geneset.xlsx")
gene_list <- toupper(trimws(as.character(geneset[[1]])))

# ============================================================
# 3. Subset to gene set genes that exist in the TPM matrix
# ============================================================
genes_in_order <- gene_list[gene_list %in% rownames(tpm_df)]
expr_subset    <- tpm_df[genes_in_order, ]

# ============================================================
# 4. Define developmental groups and assign samples
# ============================================================
groups <- c("Oocyte", "Zygote", "2-cell embryo", "4-cell embryo",
            "8-cell embryo", "Morulae", "Late blastocyst",
            "hESC passage#0", "hESC passage#10")

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
# 5. Compute per-group averages
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
# 6. log2(TPM+1) and row-scale (z-score)
# ============================================================
expr_log <- log2(expr_avg + 1)
expr_log <- t(scale(t(expr_log)))
expr_log <- expr_log[rowSums(is.na(expr_log)) == 0, ]
expr_log <- expr_log[, colSums(is.na(expr_log)) == 0]

# ============================================================
# 7. Shared color palette
# ============================================================
col_fun_vec <- colorRampPalette(c(
  "#313695", "#4575B4", "#74ADD1",
  "#FFFFFF",
  "#FDAE61", "#F46D43", "#A50026"
))(100)

# ============================================================
# 8. pheatmap output
# ============================================================
png("Embryo_Zygote-2cell_clustered.png",
    width = 1000, height = 1000, res = 180)

pheatmap(
  expr_log,
  cluster_rows  = TRUE,
  cluster_cols  = FALSE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  color         = col_fun_vec,
  main          = "Geneset Heatmap (pheatmap, vivid colors)",
  angle_col     = 90,
  cellwidth     = 16
)

dev.off()
