# Quant.sf to Gene-level TPM Matrix Pipeline
# -------------------------------------------
# Converts Salmon transcript-level quantifications into a gene-level TPM matrix
# annotated with HGNC symbols.
#

# This script:
#   1. Imports Salmon quantifications via tximport (gene-level)
#   2. Maps Ensembl gene IDs to HGNC symbols
#   3. Exports a clean gene_symbol x sample TPM matrix

# ============================================================
# 1. Install and load required packages
# ============================================================
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("tximport",        ask = FALSE)
BiocManager::install("biomaRt",         ask = FALSE)
BiocManager::install("GenomicFeatures", ask = FALSE)
install.packages("readr")
install.packages("dplyr")

library(tximport)
library(GenomicFeatures)
library(biomaRt)
library(readr)
library(dplyr)

# ============================================================
# 2. Define file paths
# ============================================================
samples <- c("Primed_1", "Primed_2", "DMSO_1", "DMSO_2",
             "PY60_24h_1", "PY60_24h_2", "PY60_48h_1", "PY60_48h_2",
             "PY60_Sorted_1", "PY60_Sorted_2", "DUX4_OE_1", "DUX4_OE_2")

files <- file.path("path/to/salmon/output", samples, "quant.sf")
names(files) <- samples

# ============================================================
# 3. Build tx2gene mapping
# ============================================================
# Either build from a GTF (preferred for reproducibility)...
txdb <- makeTxDbFromGFF("path/to/gencode.vXX.primary_assembly.annotation.gtf")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, keys = k, keytype = "TXNAME", columns = "GENEID")

# ...or load a pre-built CSV with two columns: transcript_id, gene_id
# tx2gene <- read_csv("tx2gene.csv", col_names = c("TXNAME", "GENEID"))

# ============================================================
# 4. Import gene-level TPM via tximport
# ============================================================
# By default (countsFromAbundance = "no"), tximport sums transcript-level
# TPMs into gene-level TPMs in `txi$abundance`. This is the correct
# gene-level TPM matrix.
txi <- tximport(files,
                type            = "salmon",
                tx2gene         = tx2gene,
                ignoreTxVersion = TRUE)

# Gene-level TPM matrix
gene_tpm_ensg <- as.data.frame(txi$abundance)
gene_tpm_ensg$ensembl_gene_id <- gsub("\\..*", "", rownames(gene_tpm_ensg))

# ============================================================
# 5. Map Ensembl gene IDs to HGNC symbols
# ============================================================
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

annot <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
               filters    = "ensembl_gene_id",
               values     = unique(gene_tpm_ensg$ensembl_gene_id),
               mart       = ensembl)

annot <- annot[annot$hgnc_symbol != "", ]

# Merge mapping with TPM matrix
gene_tpm <- merge(gene_tpm_ensg, annot,
                  by = "ensembl_gene_id",
                  all.x = FALSE)

# ============================================================
# 6. Collapse to one row per HGNC symbol (sum TPMs)
# ============================================================
# When multiple Ensembl genes map to the same HGNC symbol, sum their TPMs.
gene_tpm_final <- gene_tpm %>%
  group_by(hgnc_symbol, gene_biotype) %>%
  summarise(across(all_of(samples), ~ sum(.x, na.rm = TRUE)),
            .groups = "drop")

# ============================================================
# 7. Export
# ============================================================
write.table(gene_tpm_final,
            file = "TPM_gene_level_combined.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
