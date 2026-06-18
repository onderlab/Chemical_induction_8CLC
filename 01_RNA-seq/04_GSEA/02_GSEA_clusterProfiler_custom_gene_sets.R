# GSEA Analysis using clusterProfiler
# ------------------------------------
# Two workflows:
#   A) Ranking by signed -log10(p-value) using a custom gene set (Excel)
#   B) Ranking by log2FoldChange using a GMT file

library(readxl)
library(clusterProfiler)
library(enrichplot)
library(writexl)
library(ggplot2)

# ============================================================================
# WORKFLOW A: CUSTOM GENE SET FROM EXCEL + SIGNED -LOG10(P) RANKING
# ============================================================================

# ---- 1. Load DESeq2 results ----
de_res <- readxl::read_excel(
  "path/to/DESeq2_PY60_48h_vs_DMSO_results_with_HGNC_Protein_Coding_only.xlsx"
)
colnames(de_res) <- tolower(colnames(de_res))  # lowercase column names

# ---- 2. Load custom gene set (e.g. from an Excel single-column file) ----
geneset_vec <- readxl::read_excel(
  "path/to/YAP-bound_geneset_Zhao_et_al.xlsx",
  col_names = FALSE
)[[1]]
geneset_vec  <- toupper(as.character(geneset_vec))
geneset_name <- "Yap_bound_genes (Zhao et al.)"

TERM2GENE_custom <- data.frame(
  gs_name     = geneset_name,
  gene_symbol = geneset_vec
)

# ---- 3. Filter ----
de_res <- de_res[!is.na(de_res$pvalue) & !is.na(de_res$hgnc) & de_res$hgnc != "", ]
de_res$hgnc <- toupper(de_res$hgnc)
de_res <- de_res[is.finite(de_res$log2foldchange), ]
de_res <- de_res[!duplicated(de_res$hgnc), ]

# ---- 4. Build ranking metric: signed -log10(p) ----
gene_ranking <- -log10(de_res$pvalue) * sign(de_res$log2foldchange)
names(gene_ranking) <- de_res$hgnc
gene_ranking <- sort(gene_ranking, decreasing = TRUE)

# ---- 5. Run GSEA ----
gsea_res <- GSEA(
  geneList     = gene_ranking,
  TERM2GENE    = TERM2GENE_custom,
  minGSSize    = 1,
  maxGSSize    = Inf,
  pvalueCutoff = 0.25,
  nPermSimple  = 20000,
  verbose      = TRUE
)

# ---- 6. Save results ----
write_xlsx(
  gsea_res@result,
  "GSEA_YAP_Zhao_pvalueRank_results_py60_48hvsdmso.xlsx"
)

# ---- 7. Enrichment plot ----
if (nrow(gsea_res@result) > 0) {
  gs_id <- gsea_res@result$ID[1]
  p <- gseaplot2(gsea_res, geneSetID = gs_id)
  ggsave(
    "GSEA_YAP_Zhao_pvalueRank_enrichmentplot_py60_48hvsdmso.pdf",
    plot = p, width = 8, height = 6
  )
}

# ---- 8. Core enriched genes ----
if (nrow(gsea_res@result) > 0 && !is.na(gsea_res@result$core_enrichment[1])) {
  core_genes <- unique(unlist(strsplit(gsea_res@result$core_enrichment, "/")))
  cat("Core enriched genes:\n")
  print(core_genes)
  writeLines(
    core_genes,
    "GSEA_core_enriched_genes_YAP_Zhao_py60_48hvsdmso.txt"
  )
} else {
  cat("No enriched genes found.\n")
}


# ============================================================================
# WORKFLOW B: GMT FILE + LOG2FOLDCHANGE RANKING
# ============================================================================

library(readxl)
library(clusterProfiler)
library(enrichplot)
library(writexl)
library(ggplot2)

# ---- 1. Load DESeq2 results ----
de_res <- readxl::read_excel(
  "path/to/DESeq2_PY60_Sorted_vs_DMSO_results_with_HGNC_Protein_Coding_only.xlsx"
)
colnames(de_res) <- tolower(colnames(de_res))

# ---- 2. Load gene set from a GMT file ----
gmt_path <- "path/to/CORDENONSI_YAP_CONSERVED_SIGNATURE.v2025.1.Hs.gmt"
TERM2GENE_custom <- clusterProfiler::read.gmt(gmt_path)

# Optional: uppercase gene names
TERM2GENE_custom$gene <- toupper(TERM2GENE_custom$gene)

# ---- 3. Filter ----
de_res <- de_res[!is.na(de_res$log2foldchange) & !is.na(de_res$hgnc) & de_res$hgnc != "", ]
de_res$hgnc <- toupper(de_res$hgnc)
de_res <- de_res[!duplicated(de_res$hgnc), ]

# ---- 4. Build ranking metric: log2FoldChange ----
gene_ranking <- de_res$log2foldchange
names(gene_ranking) <- de_res$hgnc
gene_ranking <- sort(gene_ranking, decreasing = TRUE)

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
write_xlsx(
  gsea_res@result,
  "GSEA_YAP_Zhao_log2FC_results_py60_sortedvsdmso.xlsx"
)

# ---- 7. Enrichment plot ----
if (nrow(gsea_res@result) > 0) {
  gs_id <- gsea_res@result$ID[1]
  p <- gseaplot2(gsea_res, geneSetID = gs_id)
  ggsave(
    "GSEA_YAP_Zhao_log2FC_enrichmentplot_py60_sortedvsdmso.pdf",
    plot = p, width = 8, height = 6
  )
}

# ---- 8. Core enriched genes ----
if (nrow(gsea_res@result) > 0 && !is.na(gsea_res@result$core_enrichment[1])) {
  core_genes <- unique(unlist(strsplit(gsea_res@result$core_enrichment, "/")))
  cat("Core enriched genes:\n")
  print(core_genes)
  writeLines(
    core_genes,
    "GSEA_core_enriched_genes_YAP_Zhao_py60_sortedvsdmso.txt"
  )
} else {
  cat("No enriched genes found.\n")
}
