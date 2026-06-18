# GSEA from DESeq2 Results Using a Custom Gene Set
# --------------------------------------------------
# Two workflows:
#   A) Custom gene set provided as a single-column Excel file
#   B) Gene set(s) provided as a GMT file
#


# ============================================================
# 1. Install and load packages
# ============================================================
if (!requireNamespace("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler")
if (!requireNamespace("enrichplot",      quietly = TRUE)) BiocManager::install("enrichplot")
if (!requireNamespace("readxl",          quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("writexl",         quietly = TRUE)) install.packages("writexl")

library(clusterProfiler)
library(enrichplot)
library(readxl)
library(writexl)
library(ggplot2)
library(dplyr)


# ============================================================================
# WORKFLOW A: CUSTOM GENE SET (EXCEL) + DESeq2 STAT-BASED RANKING
# ============================================================================

# ---- 1. Load DESeq2 results ----
de_res <- readxl::read_excel(
  "path/to/DESeq2_PY60_48h_vs_DMSO_results_with_HGNC_Protein_Coding_only.xlsx"
)
colnames(de_res) <- tolower(colnames(de_res))

# ---- 2. Filter and de-duplicate ----
de_res <- de_res %>%
  filter(!is.na(stat), !is.na(hgnc), hgnc != "") %>%
  mutate(hgnc = toupper(hgnc)) %>%
  group_by(hgnc) %>%
  slice_max(abs(stat), n = 1, with_ties = FALSE) %>%
  ungroup()

# ---- 3. Build ranking from Wald statistic ----
gene_ranking <- setNames(de_res$stat, de_res$hgnc)
gene_ranking <- sort(gene_ranking, decreasing = TRUE)

# ---- 4. Load custom gene set ----
geneset_vec  <- readxl::read_excel(
  "path/to/8C-Morula_geneset_final.xlsx",
  col_names = FALSE
)[[1]]
geneset_vec  <- toupper(as.character(geneset_vec))
geneset_name <- "8C-embryo/Morula geneset"

TERM2GENE_custom <- data.frame(
  gs_name     = geneset_name,
  gene_symbol = geneset_vec
)

# ---- 5. Run GSEA ----
gsea_res <- GSEA(
  geneList     = gene_ranking,
  TERM2GENE    = TERM2GENE_custom,
  minGSSize    = 1,
  maxGSSize    = Inf,
  pvalueCutoff = 0.25,
  verbose      = FALSE
)

# ---- 6. Save results ----
write_xlsx(gsea_res@result, "GSEA_CustomSet_results.xlsx")

# ---- 7. Enrichment plot ----
if (nrow(gsea_res@result) > 0) {
  p <- gseaplot2(gsea_res, geneSetID = gsea_res@result$ID[1])
  ggsave("GSEA_CustomSet_enrichmentplot.pdf", plot = p, width = 8, height = 6)

  # ---- 8. Core enriched genes ----
  core_genes <- unique(unlist(strsplit(gsea_res@result$core_enrichment, "/")))
  print(core_genes)
} else {
  cat("No enriched gene sets at the specified cutoff.\n")
}


# ============================================================================
# WORKFLOW B: GMT FILE + DESeq2 STAT-BASED RANKING
# ============================================================================

library(clusterProfiler)
library(readxl)
library(writexl)
library(enrichplot)
library(ggplot2)
library(dplyr)

# ---- 1. Load DESeq2 results ----
de_res <- readxl::read_excel(
  "path/to/DESeq2_PY60_Sorted_vs_DMSO_results_with_HGNC_Protein_Coding_only.xlsx"
)
colnames(de_res) <- tolower(colnames(de_res))

# ---- 2. Filter and de-duplicate ----
de_res <- de_res %>%
  filter(!is.na(stat), !is.na(hgnc), hgnc != "") %>%
  mutate(hgnc = toupper(hgnc)) %>%
  group_by(hgnc) %>%
  slice_max(abs(stat), n = 1, with_ties = FALSE) %>%
  ungroup()

# ---- 3. Build ranking from Wald statistic ----
gene_ranking <- setNames(de_res$stat, de_res$hgnc)
gene_ranking <- sort(gene_ranking, decreasing = TRUE)

# ---- 4. Load gene set from a GMT file ----
gmt_path <- "path/to/YAP_TAZ_UP.v2025.1.Hs.gmt"
gmt_list <- clusterProfiler::read.gmt(gmt_path)
gmt_list$gene <- toupper(gmt_list$gene)

# ---- 5. Run GSEA ----
gsea_res <- GSEA(
  geneList     = gene_ranking,
  TERM2GENE    = gmt_list,
  minGSSize    = 10,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

# ---- 6. Save results ----
write_xlsx(gsea_res@result, "GSEA_GMT_results.xlsx")

# ---- 7. Dotplot ----
pdf("GSEA_dotplot.pdf", width = 8, height = 6)
print(dotplot(gsea_res, showCategory = 20))
dev.off()

# ---- 8. Enrichment plot per enriched gene set ----
enriched_sets <- gsea_res@result$ID[
  !is.na(gsea_res@result$NES) & gsea_res@result$p.adjust < 0.05
]

if (length(enriched_sets) > 0) {
  for (gs in enriched_sets) {
    p <- gseaplot2(gsea_res, geneSetID = gs)
    ggsave(paste0("GSEA_", gs, "_enrichmentplot.pdf"),
           plot = p, width = 8, height = 6)
  }
} else {
  cat("No enriched/plotable gene sets found.\n")
}

print("Enriched/plotted gene sets:")
print(enriched_sets)
