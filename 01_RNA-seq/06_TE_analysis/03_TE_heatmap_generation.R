# TE Heatmap Generation from TPM Values
# --------------------------------------


library(readxl)
library(dplyr)
library(pheatmap)
library(grid)

# ============================================================
# 1. Load TPM data
# ============================================================
my_data <- read_excel("path/to/Samples_TE_TPM.xlsx")

# ============================================================
# 2. Define TE list to plot
# ============================================================
te_list <- c(
  "SATR1", "SVA_D", "SVA_C", "SVA_E", "LTR14C", "SVA_B",
  "LTR5_Hs", "SVA_A", "MER34C",
  "MER9a2", "HERVK9-int", "L1_14", "MER61C", "MER11B",
  "LTR10E", "MER68B", "LTR24B", "MSTB", "MSTD-int",
  "MER11C", "MER9a1", "HERVK11-int", "HSMAR1", "LTR7B",
  "LTR8", "LTR47B", "LTR9", "LTR75_1", "MER6TE", "MST1C",
  "LTR2A2", "REP522", "MLT2A1", "MLT2A2", "HERVK13-int", "SAR"
)

# ============================================================
# 3. Filter for selected TEs
# ============================================================
my_filtered <- my_data %>% filter(Name %in% te_list)

print(paste("Number of TEs found:", nrow(my_filtered)))
print("Column names:")
print(colnames(my_filtered))

# Convert to data frame and set rownames
my_filtered <- as.data.frame(my_filtered)
rownames(my_filtered) <- my_filtered$Name
my_filtered <- my_filtered[, -1]

# ============================================================
# 4. Filter out low-expression TEs (TPM > 1 in at least one sample)
# ============================================================
my_filtered <- my_filtered[rowSums(my_filtered > 1) > 0, ]
print(paste("TEs after filtering:", nrow(my_filtered)))

# ============================================================
# 5. Reorder columns: DMSO -> PY60_48h -> PY60_Sorted
# ============================================================
column_order <- c(
  "DMSO_48h_1",    "DMSO_48h_2",
  "PY60_48h_1",    "PY60_48h_2",
  "PY60_Sorted_1", "PY60_Sorted_2"
)
available_cols <- column_order[column_order %in% colnames(my_filtered)]
my_filtered   <- my_filtered[, available_cols]

# ============================================================
# 6. Log10 transformation (standard for TPM)
# ============================================================
my_log <- log10(my_filtered + 1)

# ============================================================
# 7. Preview heatmap
# ============================================================
pheatmap(
  my_log,
  scale         = "row",
  cluster_rows  = TRUE,
  cluster_cols  = FALSE,
  color         = colorRampPalette(c(
    "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
    "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"
  ))(100),
  main          = "",
  fontsize_row  = 10,
  fontsize_col  = 11,
  fontfamily    = "Arial",
  border_color  = "grey90",
  cellwidth     = 30,
  cellheight    = 18,
  treeheight_row = 40,
  treeheight_col = 40
)

# ============================================================
# 8. Save as PDF
# ============================================================
pdf("TE_heatmap_my_samples_publication.pdf",
    width = 7, height = 12, family = "Arial")
pheatmap(
  my_log,
  scale         = "row",
  cluster_rows  = TRUE,
  cluster_cols  = FALSE,
  color         = colorRampPalette(c(
    "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
    "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"
  ))(100),
  main          = "",
  fontsize_row  = 10,
  fontsize_col  = 11,
  fontfamily    = "Arial",
  border_color  = "grey90",
  cellwidth     = 30,
  cellheight    = 18,
  treeheight_row = 40,
  treeheight_col = 40,
  legend        = TRUE
)
dev.off()

# ============================================================
# 9. Save as TIFF (300 dpi)
# ============================================================
tiff("TE_heatmap_my_samples_publication.tiff",
     width = 7, height = 12, units = "in", res = 300, family = "Arial")
pheatmap(
  my_log,
  scale         = "row",
  cluster_rows  = TRUE,
  cluster_cols  = FALSE,
  color         = colorRampPalette(c(
    "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
    "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"
  ))(100),
  main          = "",
  fontsize_row  = 10,
  fontsize_col  = 11,
  fontfamily    = "Arial",
  border_color  = "grey90",
  cellwidth     = 30,
  cellheight    = 18,
  treeheight_row = 40,
  treeheight_col = 40
)
dev.off()

print("Publication-quality heatmap generated.")
print("PDF:  TE_heatmap_my_samples_publication.pdf")
print("TIFF: TE_heatmap_my_samples_publication.tiff")
