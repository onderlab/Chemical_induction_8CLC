#!/usr/bin/env Rscript

# ============================================================================
# GENE SET ENRICHMENT ANALYSIS (GSEA)

# ============================================================================

library(fgsea)
library(msigdbr)
library(ggplot2)
library(dplyr)
library(data.table)
library(openxlsx)

cat("==============================================================\n")
cat("            GENE SET ENRICHMENT ANALYSIS (GSEA)               \n")
cat("==============================================================\n\n")

# ============================================================================
# PARAMETERS
# ============================================================================

# DESeq2 results file
deseq2_results <- "path/to/DESeq2_PY60_Sorted_vs_DMSO_all_genes.xlsx"

# Output directory
output_dir <- "path/to/GSEA_Results_PY60Sorted_vs_DMSO"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# GSEA parameters
fdr_cutoff   <- 0.25   # GSEA standard threshold
nPermSimple  <- 10000

cat("Parameters:\n")
cat(sprintf("  DESeq2 results: %s\n", deseq2_results))
cat(sprintf("  Output dir:     %s\n", output_dir))
cat(sprintf("  FDR cutoff:     %.2f\n", fdr_cutoff))
cat(sprintf("  Permutations:   %d\n\n", nPermSimple))

# ============================================================================
# LOAD DESeq2 RESULTS
# ============================================================================

cat("=== Loading DESeq2 results ===\n")

res_df <- read.xlsx(deseq2_results)
cat(sprintf("Loaded: %d genes\n", nrow(res_df)))

# Filter NA values
res_df <- res_df[!is.na(res_df$stat) & !is.na(res_df$HGNC), ]
cat(sprintf("After filtering: %d genes (NAs removed)\n\n", nrow(res_df)))

# ============================================================================
# BUILD RANKED GENE LIST
# ============================================================================

cat("=== Building ranked gene list ===\n")

gene_list <- res_df$stat
names(gene_list) <- res_df$HGNC

# Handle duplicates: keep the entry with the largest |stat| per HGNC symbol
gene_df <- data.frame(gene_name = names(gene_list), stat = gene_list)
gene_df <- gene_df %>%
  group_by(gene_name) %>%
  slice_max(abs(stat), n = 1, with_ties = FALSE) %>%
  ungroup()

gene_list <- setNames(gene_df$stat, gene_df$gene_name)
gene_list <- sort(gene_list, decreasing = TRUE)

cat(sprintf("Total genes:   %d\n", length(gene_list)))
cat(sprintf("Highest stat:  %.2f (%s)\n", gene_list[1], names(gene_list)[1]))
cat(sprintf("Lowest stat:   %.2f (%s)\n\n",
            gene_list[length(gene_list)],
            names(gene_list)[length(gene_list)]))

# Save ranked list
write.table(
  data.frame(gene = names(gene_list), stat = gene_list),
  file.path(output_dir, "ranked_gene_list.txt"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

# ============================================================================
# LOAD MSigDB GENE SETS
# ============================================================================

cat("=== Loading MSigDB gene sets ===\n")

# Hallmark
hallmark      <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_list <- split(hallmark$gene_symbol, hallmark$gs_name)
cat(sprintf("Hallmark: %d gene sets\n", length(hallmark_list)))

# Reactome
reactome      <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
reactome_list <- split(reactome$gene_symbol, reactome$gs_name)
cat(sprintf("Reactome: %d pathways\n", length(reactome_list)))

# GO Biological Process
go_bp      <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
go_bp_list <- split(go_bp$gene_symbol, go_bp$gs_name)
cat(sprintf("GO BP:    %d terms\n\n", length(go_bp_list)))

# ============================================================================
# fgsea - HALLMARK
# ============================================================================

cat("=== GSEA: Hallmark Gene Sets ===\n")

set.seed(42)
fgsea_hallmark <- fgsea(
  pathways    = hallmark_list,
  stats       = gene_list,
  minSize     = 15,
  maxSize     = 500,
  nPermSimple = nPermSimple
)
fgsea_hallmark <- fgsea_hallmark[order(fgsea_hallmark$pval), ]

cat(sprintf("Total: %d pathways\n", nrow(fgsea_hallmark)))
cat(sprintf("Significant (padj < %.2f): %d\n", fdr_cutoff,
            sum(fgsea_hallmark$padj < fdr_cutoff, na.rm = TRUE)))

hallmark_up   <- fgsea_hallmark[fgsea_hallmark$NES > 0 & fgsea_hallmark$padj < fdr_cutoff, ]
hallmark_down <- fgsea_hallmark[fgsea_hallmark$NES < 0 & fgsea_hallmark$padj < fdr_cutoff, ]

cat(sprintf("  UP in PY60:   %d\n", nrow(hallmark_up)))
cat(sprintf("  DOWN in PY60: %d\n", nrow(hallmark_down)))

# ============================================================================
# fgsea - REACTOME
# ============================================================================

cat("\n=== GSEA: Reactome Pathways ===\n")

set.seed(42)
fgsea_reactome <- fgsea(
  pathways    = reactome_list,
  stats       = gene_list,
  minSize     = 15,
  maxSize     = 500,
  nPermSimple = nPermSimple
)
fgsea_reactome <- fgsea_reactome[order(fgsea_reactome$pval), ]

cat(sprintf("Significant (padj < %.2f): %d / %d\n",
            fdr_cutoff,
            sum(fgsea_reactome$padj < fdr_cutoff, na.rm = TRUE),
            nrow(fgsea_reactome)))

# ============================================================================
# fgsea - GO BP (optional, slow)
# ============================================================================

cat("\n=== GSEA: GO Biological Process ===\n")
cat("(This may take several minutes...)\n")

set.seed(42)
fgsea_go_bp <- fgsea(
  pathways    = go_bp_list,
  stats       = gene_list,
  minSize     = 15,
  maxSize     = 500,
  nPermSimple = nPermSimple
)
fgsea_go_bp <- fgsea_go_bp[order(fgsea_go_bp$pval), ]

cat(sprintf("Significant (padj < %.2f): %d / %d\n",
            fdr_cutoff,
            sum(fgsea_go_bp$padj < fdr_cutoff, na.rm = TRUE),
            nrow(fgsea_go_bp)))

# ============================================================================
# PRINT TOP RESULTS
# ============================================================================

cat("\n==============================================================\n")
cat("                    TOP ENRICHED PATHWAYS                     \n")
cat("==============================================================\n")

cat("\n=== TOP 10 HALLMARK (UP in PY60) ===\n")
if (nrow(hallmark_up) > 0) {
  top_up <- head(hallmark_up[order(hallmark_up$NES, decreasing = TRUE), ], 10)
  for (i in 1:nrow(top_up)) {
    name <- gsub("HALLMARK_", "", top_up$pathway[i])
    cat(sprintf("%2d. %-40s NES=%5.2f  padj=%.3f\n",
                i, name, top_up$NES[i], top_up$padj[i]))
  }
} else {
  cat("  None\n")
}

cat("\n=== TOP 10 HALLMARK (DOWN in PY60) ===\n")
if (nrow(hallmark_down) > 0) {
  top_down <- head(hallmark_down[order(hallmark_down$NES), ], 10)
  for (i in 1:nrow(top_down)) {
    name <- gsub("HALLMARK_", "", top_down$pathway[i])
    cat(sprintf("%2d. %-40s NES=%5.2f  padj=%.3f\n",
                i, name, top_down$NES[i], top_down$padj[i]))
  }
} else {
  cat("  None\n")
}

cat("\n=== TOP 15 REACTOME PATHWAYS ===\n")
reactome_sig <- fgsea_reactome[fgsea_reactome$padj < fdr_cutoff, ]
if (nrow(reactome_sig) > 0) {
  top_reactome <- head(reactome_sig[order(abs(reactome_sig$NES), decreasing = TRUE), ], 15)
  for (i in 1:nrow(top_reactome)) {
    name <- gsub("REACTOME_", "", top_reactome$pathway[i])
    name <- substr(name, 1, 50)
    dir  <- ifelse(top_reactome$NES[i] > 0, "UP", "DOWN")
    cat(sprintf("%2d. %-50s %s NES=%5.2f padj=%.3f\n",
                i, name, dir, top_reactome$NES[i], top_reactome$padj[i]))
  }
} else {
  cat("  No significant pathways found\n")
}

# ============================================================================
# VISUALIZATION
# ============================================================================

cat("\n=== Generating plots ===\n")

# Enrichment plots (top 6 Hallmark)
pdf(file.path(output_dir, "enrichment_plots.pdf"), width = 10, height = 6)
if (nrow(hallmark_up) > 0) {
  for (pathway in head(hallmark_up$pathway, 6)) {
    p <- plotEnrichment(hallmark_list[[pathway]], gene_list) +
      labs(title = gsub("HALLMARK_", "", pathway)) +
      theme_minimal()
    print(p)
  }
}
if (nrow(hallmark_down) > 0) {
  for (pathway in head(hallmark_down$pathway, 6)) {
    p <- plotEnrichment(hallmark_list[[pathway]], gene_list) +
      labs(title = gsub("HALLMARK_", "", pathway)) +
      theme_minimal()
    print(p)
  }
}
dev.off()
cat("Enrichment plots written\n")

# Barplot - Hallmark
hallmark_sig <- fgsea_hallmark[fgsea_hallmark$padj < fdr_cutoff, ]
if (nrow(hallmark_sig) > 0) {
  hallmark_sig$pathway_clean <- gsub("HALLMARK_", "", hallmark_sig$pathway)
  hallmark_sig$pathway_clean <- gsub("_", " ", hallmark_sig$pathway_clean)
  top20 <- head(hallmark_sig[order(abs(hallmark_sig$NES), decreasing = TRUE), ], 20)
  top20$pathway_clean <- factor(top20$pathway_clean,
                                levels = top20$pathway_clean[order(top20$NES)])

  p_bar <- ggplot(top20, aes(x = NES, y = pathway_clean, fill = NES > 0)) +
    geom_col() +
    scale_fill_manual(
      values = c("TRUE" = "#B2182B", "FALSE" = "#2166AC"),
      labels = c("TRUE" = "UP in PY60", "FALSE" = "DOWN in PY60"),
      name   = ""
    ) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(
      title    = "Top 20 Hallmark Pathways",
      subtitle = sprintf("GSEA: PY60_Sorted vs DMSO (FDR < %.2f)", fdr_cutoff),
      x = "Normalized Enrichment Score (NES)",
      y = ""
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle   = element_text(hjust = 0.5),
      legend.position = "bottom"
    )
  ggsave(file.path(output_dir, "hallmark_barplot.png"),
         p_bar, width = 10, height = 8, dpi = 300)
  cat("Barplot written\n")

  # Dotplot
  top15 <- head(hallmark_sig[order(hallmark_sig$pval), ], 15)
  top15$pathway_clean <- factor(top15$pathway_clean,
                                levels = top15$pathway_clean[order(top15$NES)])
  p_dot <- ggplot(top15, aes(x = NES, y = pathway_clean)) +
    geom_point(aes(size = -log10(padj), color = NES)) +
    scale_color_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, name = "NES"
    ) +
    scale_size_continuous(name = "-log10(FDR)") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(
      title    = "Top 15 Hallmark Pathways",
      subtitle = "Sized by significance, colored by direction",
      x = "Normalized Enrichment Score",
      y = ""
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )
  ggsave(file.path(output_dir, "hallmark_dotplot.png"),
         p_dot, width = 10, height = 7, dpi = 300)
  cat("Dotplot written\n")
}

# ============================================================================
# SAVE RESULTS
# ============================================================================

cat("\n=== Saving results ===\n")

# Hallmark
write.table(
  fgsea_hallmark[, c("pathway", "pval", "padj", "ES", "NES", "size")],
  file.path(output_dir, "GSEA_hallmark_all.txt"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
if (nrow(hallmark_sig) > 0) {
  write.table(
    hallmark_sig[, c("pathway", "pval", "padj", "ES", "NES", "size")],
    file.path(output_dir, "GSEA_hallmark_significant.txt"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
}

# Reactome
write.table(
  fgsea_reactome[, c("pathway", "pval", "padj", "ES", "NES", "size")],
  file.path(output_dir, "GSEA_reactome_all.txt"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
if (nrow(reactome_sig) > 0) {
  write.table(
    reactome_sig[, c("pathway", "pval", "padj", "ES", "NES", "size")],
    file.path(output_dir, "GSEA_reactome_significant.txt"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
}

# GO BP
go_bp_sig <- fgsea_go_bp[fgsea_go_bp$padj < fdr_cutoff, ]
write.table(
  fgsea_go_bp[, c("pathway", "pval", "padj", "ES", "NES", "size")],
  file.path(output_dir, "GSEA_GO_BP_all.txt"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
if (nrow(go_bp_sig) > 0) {
  write.table(
    go_bp_sig[, c("pathway", "pval", "padj", "ES", "NES", "size")],
    file.path(output_dir, "GSEA_GO_BP_significant.txt"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
}

cat("All results saved\n")

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n==============================================================\n")
cat("                    GSEA COMPLETE                             \n")
cat("==============================================================\n\n")

cat("RESULTS:\n")
cat(sprintf("  Hallmark:  %d / %d significant\n",
            sum(fgsea_hallmark$padj < fdr_cutoff, na.rm = TRUE),
            nrow(fgsea_hallmark)))
cat(sprintf("  Reactome:  %d / %d significant\n",
            sum(fgsea_reactome$padj < fdr_cutoff, na.rm = TRUE),
            nrow(fgsea_reactome)))
cat(sprintf("  GO BP:     %d / %d significant\n",
            sum(fgsea_go_bp$padj < fdr_cutoff, na.rm = TRUE),
            nrow(fgsea_go_bp)))

cat("\nOUTPUT:\n")
cat(sprintf("  Output directory: %s\n", output_dir))
cat("\nFiles:\n")
cat("  - hallmark_barplot.png\n")
cat("  - hallmark_dotplot.png\n")
cat("  - enrichment_plots.pdf\n")
cat("  - GSEA_hallmark_significant.txt\n")
cat("  - GSEA_reactome_significant.txt\n")
cat("  - GSEA_GO_BP_significant.txt\n")
cat("  - ranked_gene_list.txt\n")

cat("\nINTERPRETATION:\n")
cat("  NES > 0  = Pathway UP-regulated in PY60\n")
cat("  NES < 0  = Pathway DOWN-regulated in PY60\n")

cat("\n==============================================================\n\n")
