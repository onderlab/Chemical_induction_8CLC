# GSEA Plot: NES vs FDR (q-value) with highlighted gene set
# ----------------------------------------------------------


library(readxl)
library(dplyr)
library(clusterProfiler)
library(GSEABase)
library(ggplot2)
library(openxlsx)

# ============================================================
# 1. Load DESeq2 results
# ============================================================
de_res <- readxl::read_excel(
  "path/to/DESeq2_PY60_48h_vs_DMSO_results_with_HGNC_Protein_Coding_only.xlsx"
)
colnames(de_res) <- tolower(colnames(de_res))

# ============================================================
# 2. Build GSEA ranking metric: sign(log2FC) * -log10(padj)
# ============================================================
# Adjust column name if needed: 'hgnc', 'hgnc_symbol', 'gene_symbol', etc.
symbol_col <- "hgnc"
stopifnot(symbol_col %in% colnames(de_res))

gsea_input <- de_res %>%
  filter(!is.na(!!sym(symbol_col)),
         !!sym(symbol_col) != "",
         !is.na(padj)) %>%
  mutate(ranking_score = sign(log2foldchange) *
                         -log10(ifelse(padj == 0, 1e-300, padj))) %>%
  filter(is.finite(ranking_score), !duplicated(!!sym(symbol_col)))

gene_ranking <- gsea_input$ranking_score
names(gene_ranking) <- toupper(gsea_input[[symbol_col]])
gene_ranking <- sort(gene_ranking, decreasing = TRUE)

# ============================================================
# 3. Load all GMT files from a folder
# ============================================================
gmt_folder <- "path/to/msigdb_gmt"
gmt_files  <- list.files(gmt_folder, pattern = "\\.gmt$", full.names = TRUE)

gmt_list      <- lapply(gmt_files, getGmt)
all_gmt_flat  <- unlist(gmt_list, recursive = FALSE)

gmt_df <- do.call(rbind, lapply(all_gmt_flat, function(gs) {
  data.frame(gs_name     = gs@setName,
             gene_symbol = gs@geneIds,
             stringsAsFactors = FALSE)
}))

# ============================================================
# 4. Run GSEA on all collections
# ============================================================
gsea_res <- GSEA(
  geneList     = gene_ranking,
  TERM2GENE    = gmt_df,
  minGSSize    = 15,
  maxGSSize    = 2000,
  pvalueCutoff = 1,
  nPermSimple  = 20000
)

# ============================================================
# 5. Prepare results table
# ============================================================
gsea_df <- as.data.frame(gsea_res@result)
gsea_df <- gsea_df %>% filter(qvalue < 1)
gsea_df$qvalue[gsea_df$qvalue == 0] <- 1e-10

# ============================================================
# 6. Define a gene set to highlight
# ============================================================
custom_label <- "8C-Morula"  # must match the gene set name in the GMT
custom_point <- gsea_df %>% filter(ID == custom_label)

# ============================================================
# 7. NES vs FDR plot
# ============================================================
p <- ggplot(gsea_df, aes(x = qvalue, y = NES)) +
  geom_point(
    shape = 21, size = 1.6, stroke = 0.4,
    color = "black", fill = "skyblue2"
  ) +
  geom_segment(
    data = custom_point,
    aes(x = qvalue, xend = qvalue + 0.10, y = NES, yend = NES + 0.13),
    color = "red", linewidth = 0.7
  ) +
  geom_text(
    data = custom_point,
    aes(x = qvalue + 0.13, y = NES + 0.14, label = custom_label),
    color = "red", fontface = "bold", size = 4.2, hjust = 0
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(-3, 3)) +
  theme_bw(base_size = 13) +
  labs(
    x     = "FDR q-value (0-1, linear scale)",
    y     = "NES (Normalized Enrichment Score)",
    title = "GSEA: NES vs FDR q-value (All MSigDB + 8C-Morula gene set)"
  )

print(p)

# ============================================================
# 8. Save outputs
# ============================================================
ggsave("PY60_48h_vs_DMSO_GSEA_NES_vs_FDR_with8Cmorula.png",
       plot = p, width = 8, height = 10, dpi = 350, units = "in")

write.xlsx(gsea_df,
           file = "PY60_48h_vs_DMSO_GSEA_NES_vs_FDR_with8Cmorula.xlsx",
           rowNames = FALSE)
