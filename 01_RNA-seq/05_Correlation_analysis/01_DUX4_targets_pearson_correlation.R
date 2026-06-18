# DUX4 Target Gene Response: PY60 vs DUX4_OE Correlation
# -------------------------------------------------------

library(readxl)
library(dplyr)
library(ggplot2)
library(ggrepel)

# ============================================================
# 1. Load files
# ============================================================
dux4_targets <- read_xlsx("path/to/DUX4_targetgenes_resnicketal.xlsx")
dux4_deseq2  <- read_xlsx(
  "path/to/DESeq2_DUX4_OE_vs_DMSO_results_with_HGNC_Protein_Coding_only.xlsx"
)
py60_deseq2  <- read_xlsx(
  "path/to/DESeq2_PY60_Sorted_vs_DMSO_results_with_HGNC_Protein_Coding_only.xlsx"
)

# ============================================================
# 2. Prepare data
# ============================================================
colnames(dux4_targets)[1] <- "HGNC"
dux4_genes <- toupper(trimws(dux4_targets$HGNC))

dux4_deseq2 <- dux4_deseq2 %>% rename(HGNC = HGNC, log2fc_dux4 = log2FoldChange)
py60_deseq2 <- py60_deseq2 %>% rename(HGNC = HGNC, log2fc_py60 = log2FoldChange)

dux4_deseq2$HGNC <- toupper(dux4_deseq2$HGNC)
py60_deseq2$HGNC <- toupper(py60_deseq2$HGNC)

# ============================================================
# 3. Merge tables on DUX4 target genes
# ============================================================
merged <- inner_join(
  dux4_deseq2 %>% filter(HGNC %in% dux4_genes) %>% select(HGNC, log2fc_dux4),
  py60_deseq2 %>% filter(HGNC %in% dux4_genes) %>% select(HGNC, log2fc_py60),
  by = "HGNC"
)

# Genes to highlight
highlight_genes <- c("LEUTX", "TPRX1", "ZSCAN4", "KLF17",
                     "DUXA", "DUXB", "ZSCAN5B", "RFPL2", "SLC34A2")

merged <- merged %>%
  mutate(highlight = ifelse(HGNC %in% highlight_genes, "yes", "no"))

# ============================================================
# 4. Pearson correlation
# ============================================================
cor_test <- cor.test(merged$log2fc_dux4, merged$log2fc_py60)

# ============================================================
# 5. Custom axis breaks and labels
# ============================================================
x_breaks <- seq(-5, ceiling(max(merged$log2fc_dux4)), by = 5)
y_breaks <- seq(-5, ceiling(max(merged$log2fc_py60)), by = 5)
x_labels <- ifelse(x_breaks == -5, "", x_breaks)
y_labels <- ifelse(y_breaks == -5, "", y_breaks)

# ============================================================
# 6. Build the plot
# ============================================================
ggplot(merged, aes(x = log2fc_dux4, y = log2fc_py60)) +
  geom_point(color = "gray60", size = 2, alpha = 0.9) +
  geom_smooth(method = "lm", se = TRUE, color = "black",
              fill = "lightblue", alpha = 0.3) +

  # Dashed lines at x=0 and y=0
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +

  # Highlighted genes
  geom_point(data = merged %>% filter(highlight == "yes"),
             aes(x = log2fc_dux4, y = log2fc_py60),
             color = "blue", size = 2) +
  geom_text_repel(
    data           = merged %>% filter(highlight == "yes"),
    aes(label = HGNC),
    color          = "blue",
    size           = 3,
    fontface       = "bold",
    segment.color  = "blue",
    segment.size   = 0.6,
    segment.alpha  = 0.7
  ) +

  # Axes and theme
  coord_cartesian(xlim = c(-5, NA), ylim = c(-5, NA)) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels) +
  scale_y_continuous(breaks = y_breaks, labels = y_labels) +

  theme_minimal(base_size = 14) +
  labs(
    title    = "DUX4 Target Gene Response",
    subtitle = paste("Pearson r =", round(cor_test$estimate, 2),
                     "| p =", signif(cor_test$p.value, 2)),
    x = "log2FC (DUX4_OE vs DMSO)",
    y = "log2FC (PY60 vs DMSO)"
  ) +
  theme(
    panel.grid         = element_blank(),
    plot.title         = element_text(face = "bold"),
    plot.subtitle      = element_text(size = 12),
    axis.line          = element_line(color = "black", linewidth = 0.6),
    axis.ticks         = element_line(color = "black", linewidth = 0.4),
    axis.ticks.length  = unit(0.2, "cm")
  )

ggsave("dux4_py60_correlation_plot.pdf", width = 7, height = 7)
