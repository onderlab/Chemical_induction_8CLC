# Heatmap Generation Using VST Normalization (from tximport)
# -----------------------------------------------------------
# Pipeline:
#   1. Import Salmon quantifications via tximport
#      (countsFromAbundance = "no" -- correct setting for DESeq2)
#   2. Build DESeq2 dataset and apply VST normalization
#   3. Map Ensembl IDs to HGNC symbols
#   4. Subset to a user-defined gene set
#   5. Plot a row-scaled heatmap

library(tximport)
library(GenomicFeatures)
library(DESeq2)
library(pheatmap)
library(readxl)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)

# ============================================================
# 1. Sample names and file paths
# ============================================================
samples <- c("Primed_1", "Primed_2",
             "DMSO_1", "DMSO_2",
             "VTP_48h_1", "VTP_48h_2",
             "PY60_24h_1", "PY60_24h_2",
             "PY60_48h_1", "PY60_48h_2",
             "PY60_Sorted_1", "PY60_Sorted_2",
             "DUX4_OE_1", "DUX4_OE_2")

files <- file.path("path/to/salmon_output", samples, "quant.sf")
names(files) <- samples

# ============================================================
# 2. Build tx2gene from a GTF file
# ============================================================
txdb <- makeTxDbFromGFF("path/to/gencode.v48.primary_assembly.annotation.gtf.gz")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, keys = k, keytype = "TXNAME", columns = "GENEID")

# ============================================================
# 3. tximport: gene-level counts (DESeq2-ready)
# ============================================================
# IMPORTANT: countsFromAbundance = "no" (default) is the correct setting
# for DESeq2 input. With this option, tximport supplies:
#   - raw counts (txi$counts)
#   - an avgTxLength matrix used by DESeq2 as a length offset
# Using "lengthScaledTPM" or "scaledTPM" instead double-corrects for
# transcript length and breaks DESeq2's dispersion estimation.
txi <- tximport(files,
                type    = "salmon",
                tx2gene = tx2gene,
                countsFromAbundance = "no")

# ============================================================
# 4. Define conditions
# ============================================================
conditions <- c("Primed", "Primed",
                "DMSO", "DMSO",
                "VTP_48h", "VTP_48h",
                "PY60_24h", "PY60_24h",
                "PY60_48h", "PY60_48h",
                "PY60_Sorted", "PY60_Sorted",
                "DUX4_OE", "DUX4_OE")

coldata <- data.frame(row.names = samples,
                      condition = factor(conditions))

# ============================================================
# 5. DESeq2 dataset and group-size-aware prefilter
# ============================================================
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ condition)

# Keep genes with >=10 counts in at least as many samples as the smallest group
smallestGroupSize <- min(table(coldata$condition))
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds  <- dds[keep, ]

# ============================================================
# 6. VST normalization
# ============================================================
vsd     <- vst(dds, blind = TRUE)
vsd_mat <- assay(vsd)

# ============================================================
# 7. Map Ensembl IDs to HGNC symbols (mart_export.txt from BioMart)
# ============================================================
mapping <- read.table(
  "path/to/mart_export.txt",
  header           = TRUE,
  sep              = "\t",
  quote            = "",
  stringsAsFactors = FALSE,
  check.names      = FALSE
)
names(mapping) <- trimws(names(mapping))
names(mapping) <- gsub("\\s+", " ", names(mapping))

# Strip version suffix from Ensembl IDs
mapping$ENSG_nover <- gsub("\\..*", "", mapping$`Gene stable ID`)

# Strip version suffix from rownames of the VST matrix
ens_ids <- gsub("\\..*", "", rownames(vsd_mat))

# Build mapping
hgnc_map   <- setNames(mapping$`HGNC symbol`, mapping$ENSG_nover)
hgnc_names <- hgnc_map[ens_ids]

# Keep only mapped, non-empty symbols
vsd_mat_hgnc <- vsd_mat[!is.na(hgnc_names) & hgnc_names != "", ]
rownames(vsd_mat_hgnc) <- hgnc_names[!is.na(hgnc_names) & hgnc_names != ""]

# ============================================================
# 8. Load a gene set (single column of HGNC symbols)
# ============================================================
geneset     <- read_xlsx("path/to/8C-Morula_geneset_final.xlsx")
geneset_ids <- toupper(trimws(geneset[[1]]))

rownames(vsd_mat_hgnc) <- toupper(trimws(rownames(vsd_mat_hgnc)))

# ============================================================
# 9. Subset to the gene set
# ============================================================
vsd_sub <- vsd_mat_hgnc[rownames(vsd_mat_hgnc) %in% geneset_ids, ]

# ============================================================
# 10. Heatmap (no labels, no clustering, vivid color palette)
# ============================================================
# NOTE: VST already stabilises variance, and `scale = "row"` then computes
# z-scores PER ROW. The two combined produce strongly normalized values
# that emphasise relative (not absolute) expression across samples.
# Color intensity should be interpreted as a within-gene comparison.
png("geneset_heatmap_zscores.png",
    width = 1000, height = 1000, res = 180)

col_fun <- colorRampPalette(c(
  "#313695", "#4575B4", "#74ADD1",
  "#FFFFFF",
  "#FDAE61", "#F46D43", "#A50026"
))(100)

pheatmap(
  vsd_sub,
  cluster_rows    = FALSE,
  cluster_cols    = FALSE,
  show_rownames   = FALSE,
  show_colnames   = TRUE,
  scale           = "row",
  annotation_col  = data.frame(condition = conditions,
                               row.names = samples),
  color           = col_fun,
  main            = "Geneset Heatmap (row z-scores of VST values)"
)

dev.off()
