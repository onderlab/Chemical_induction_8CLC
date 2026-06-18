# Volcano Plot from DESeq2 Results (Protein-Coding Genes)
# --------------------------------------------------------


library(openxlsx)
library(EnhancedVolcano)
library(ggplot2)

# ============================================================
# 1. Load DESeq2 results
# ============================================================
results_gene <- read.xlsx("path/to/DESeq2_PY60_Sorted_vs_DMSO_results_with_HGNC.xlsx")

# Keep only protein-coding genes
results_gene <- subset(results_gene, gene_biotype == "protein_coding")

# ============================================================
# 2. Define marker genes to highlight
# ============================================================
marker_genes <- c(
  "ZSCAN5B", "KDM4E", "TRIM43B", "ZSCAN4", "MBD3L2", "SCL34A2",
  "TPRX1", "DUXA", "DUXB", "RFPL2", "LEUTX", "RFPL4A", "CCNA1", "KLF17"
)

# ============================================================
# 3. Get top 5 upregulated genes (by log2FoldChange among significant hits)
# ============================================================
top5_ENSG <- results_gene[results_gene$log2FoldChange > 0 & results_gene$padj < 0.05, ]
top5_ENSG <- top5_ENSG[order(top5_ENSG$log2FoldChange, decreasing = TRUE), ]
top5_ENSG <- unique(top5_ENSG$ENSG)[1:5]
top5_ENSG <- na.omit(top5_ENSG)

# Map marker gene names to ENSG IDs
marker_ENSG <- unique(results_gene$ENSG[results_gene$HGNC %in% marker_genes])
marker_ENSG <- na.omit(marker_ENSG)

# Combine label set
label_ENSG <- unique(c(top5_ENSG, marker_ENSG))

# ============================================================
# 4. Add labels to the data frame
# ============================================================
results_gene$mylabel <- ifelse(
  results_gene$ENSG %in% label_ENSG &
    !is.na(results_gene$HGNC) &
    results_gene$HGNC != "",
  results_gene$HGNC,
  ""
)

selectLab_labels <- unique(results_gene$mylabel[results_gene$mylabel != ""])
cat("Labels to be displayed:", selectLab_labels, "\n")

# ============================================================
# 5. Build the volcano plot
# ============================================================
p <- EnhancedVolcano(
  results_gene,
  lab       = results_gene$mylabel,
  selectLab = selectLab_labels,
  x         = "log2FoldChange",
  y         = "padj",
  pCutoff   = 0.05,
  FCcutoff  = 1,
  pointSize = 1,
  labSize   = 3.5,
  drawConnectors  = TRUE,
  widthConnectors = 0.5,
  colConnectors   = "black",
  arrowheads      = FALSE,
  title           = "PY60_Sorted vs DMSO (protein_coding)",
  legendLabels    = c("NS", "Log2FC", "FDR", "FDR & Log2FC"),
  legendPosition  = "right",
  col             = c("grey30", "grey30", "royalblue", "red2"),
  ylim            = c(0, 300)
)

# ============================================================
# 6. Export plot
# ============================================================
png("volcano_marker_PY60_Sorted.png",
    width = 1600, height = 1200, res = 180)
print(p + ylab("-Log10(FDR)"))
dev.off()
