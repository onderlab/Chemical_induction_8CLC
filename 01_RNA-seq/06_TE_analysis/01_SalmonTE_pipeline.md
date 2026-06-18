# SalmonTE: Transposable Element Quantification and Differential Expression

## 1. Build the SalmonTE index

```bash
path/to/SalmonTE/salmon/darwin/bin/salmon index \
  -t path/to/SalmonTE/reference/hs/hs_origin.fa \
  -i path/to/SalmonTE/reference/hs \
  -p 8
```

## 2. Run SalmonTE quantification

```bash
path/to/SalmonTE/SalmonTE.py quant \
  --reference=hs \
  --outpath path/to/TE_output \
  --num_threads 8 \
  path/to/Trimmed_Fastq_Files/
```

## 3. Subset condition table for the comparison of interest

```bash
mkdir -p TE_in_PY60_Sorted_vs_DMSO_48h

python3 - <<'PY'
import pandas as pd

expr = pd.read_csv("EXPR.csv")
cond = pd.read_csv("condition.csv", dtype=str)

# Keep only the conditions of interest (must match the labels exactly)
keep = {"DMSO_48h", "PY60_Sorted"}

cond_f = cond[cond["condition"].isin(keep)].copy()
cols = ["TE"] + cond_f["SampleID"].tolist()
expr_f = expr.loc[:, cols]

expr_f.to_csv("TE_in_PY60_Sorted_vs_DMSO_48h/EXPR.csv", index=False)
cond_f.to_csv("TE_in_PY60_Sorted_vs_DMSO_48h/condition.csv", index=False)

print("Kept samples:", len(cond_f))
print(sorted(cond_f["SampleID"].tolist()))
PY

# Quick check
head -n 2 TE_in_PY60_Sorted_vs_DMSO_48h/EXPR.csv
column -s, -t TE_in_PY60_Sorted_vs_DMSO_48h/condition.csv
```

## 4. Run SalmonTE differential expression

```bash
path/to/SalmonTE/SalmonTE.py test \
  --inpath=path/to/TE/TE_in_PY60_Sorted_vs_DMSO_48h \
  --outpath=path/to/TE/DE_PY60_Sorted_vs_DMSO_48h \
  --tabletype=xls \
  --figtype=pdf \
  --analysis_type=DE \
  --conditions="DMSO_48h,PY60_Sorted"
```

## 5. Element-level DE + volcano plot (R)

```r
# ==== PY60_Sorted vs DMSO_48h | Element-level DE + Volcano ====

suppressPackageStartupMessages({
  library(tximport); library(DESeq2); library(dplyr); library(readr)
  library(EnhancedVolcano); library(ggplot2); library(scales); library(tibble)
  library(openxlsx)
})

# ---- Paths & comparison ----
cmp_inpath <- "path/to/TE"                            # contains condition.csv
base_dir   <- "path/to/TE"                            # each sample has quant.sf here
ref_dir    <- "path/to/SalmonTE/reference/hs"         # contains clades.csv
outdir     <- "path/to/TE/DE_PY60_Sorted_vs_DMSO_48h"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

condA <- "PY60_Sorted"   # numerator
condB <- "DMSO_48h"      # denominator

# ---- Read sample metadata ----
cond  <- read.csv(file.path(cmp_inpath, "condition.csv"), stringsAsFactors = TRUE)
files <- file.path(base_dir, cond$SampleID, "quant.sf")
names(files) <- cond$SampleID

# ---- Build tx2element on the fly from clades.csv (strips _dupNNN) ----
clades <- read.csv(file.path(ref_dir, "clades.csv"))
elem_base <- function(x) sub("_dup[0-9]+$", "", as.character(x))
tx2element <- clades |> transmute(TXNAME = name, GENEID = elem_base(name)) |> distinct()
write.csv(tx2element, file.path(outdir, "tx2element.csv"), row.names = FALSE)

# ---- tximport + DESeq2 (element-level aggregation) ----
# countsFromAbundance = "no" (default) -- correct setting for DESeq2.
txi <- tximport(files,
                type    = "salmon",
                tx2gene = tx2element,
                ignoreTxVersion = TRUE,
                countsFromAbundance = "no")

dds <- DESeqDataSetFromTximport(txi, colData = cond, design = ~ condition)

# ---- Filter: drop low-count TEs (>=10 in at least n_smallestGroup samples) ----
smallestGroupSize <- min(table(cond$condition))
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds  <- dds[keep, ]

dds      <- DESeq(dds, fitType = "local")
res      <- results(dds, contrast = c("condition", condA, condB))
res_elem <- as.data.frame(res) |> rownames_to_column("element")

# ---- Keep only necessary columns ----
res_elem <- res_elem |>
  select(element, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)

# ---- Save results ----
write.csv(res_elem, file.path(outdir, "DE_element_level.csv"), row.names = FALSE)
write.xlsx(res_elem, file.path(outdir, "DE_element_level.xlsx"), rowNames = FALSE)

# ---- Volcano plot: label only significant points ----
alpha_p <- 0.05
lfc_cut <- 1
res_v <- res_elem |>
  mutate(
    pvalue_plot = ifelse(is.na(pvalue), 1, pvalue),
    label       = ifelse(!is.na(padj) & padj < alpha_p & abs(log2FoldChange) > lfc_cut,
                         element, NA)
  )

p <- EnhancedVolcano(
  res_v,
  lab       = res_v$label,
  x         = "log2FoldChange",
  y         = "pvalue_plot",
  pCutoff   = alpha_p,
  FCcutoff  = lfc_cut,
  pointSize = 2,
  labSize   = 3.5,
  drawConnectors  = TRUE,
  widthConnectors = 0.5,
  colConnectors   = "black",
  arrowheads      = FALSE,
  title           = sprintf("%s vs %s (element level)", condA, condB),
  legendLabels    = c("NS", "Log2FC", "p-value", "p & Log2FC"),
  legendPosition  = "right",
  col             = c("grey30", "grey30", "royalblue", "red2")
) +
  labs(y = expression(-log[10] ~ "(p-value)")) +
  coord_cartesian(xlim = c(-8, 8)) +
  scale_x_continuous(breaks = seq(-8, 8, 2),
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(breaks = pretty_breaks(6),
                     expand = expansion(mult = c(0.02, 0.08))) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_line(linewidth = 0.3, colour = "grey86"),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 14),
    legend.title     = element_blank()
  )

ggsave(file.path(outdir, "volcano_element.png"), p, width = 9, height = 7, dpi = 300)
ggsave(file.path(outdir, "volcano_element.pdf"), p, width = 9, height = 7)
```
