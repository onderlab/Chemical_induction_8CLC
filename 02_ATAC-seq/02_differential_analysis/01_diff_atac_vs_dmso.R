#!/usr/bin/env Rscript
# =========================================================
# diff_atac_vs_dmso.R
# ---------------------------------------------------------
# For each non-DMSO condition, runs DESeq2 vs DMSO using the
# hkTSS-derived size factors, identifies gained and lost ATAC
# peaks, and writes BED files for downstream plotHeatmap.
#
# Usage:
#   Rscript diff_atac_vs_dmso.R \
#     <consensus_peak_counts.txt> \
#     <samples.tsv> \
#     <hkTSS_scale_factors.tsv> \
#     <out_dir>
#
# Outputs (per comparison, e.g. PY60_vs_DMSO):
#   <out_dir>/PY60_vs_DMSO_gained.bed
#   <out_dir>/PY60_vs_DMSO_lost.bed
#   <out_dir>/PY60_vs_DMSO_DESeq2_results.csv
# =========================================================

args        <- commandArgs(trailingOnly = TRUE)
count_file  <- args[1]   # featureCounts output
sample_file <- args[2]   # samples.tsv
hk_sf_file  <- args[3]   # hkTSS_scale_factors.tsv
out_dir     <- args[4]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
})

# --- 1. Load data -----------------------------------------------------------
samples <- as.data.frame(fread(sample_file))
rownames(samples) <- samples$sample

fc <- as.data.frame(fread(count_file, skip = 1))
count_cols <- grep("\\.shifted\\.bam$", colnames(fc), value = TRUE)
count_mat  <- fc[, count_cols, drop = FALSE]
rownames(count_mat) <- fc$Geneid
colnames(count_mat) <- sub("\\.shifted\\.bam$", "", basename(count_cols))
count_mat  <- as.matrix(count_mat[, samples$sample, drop = FALSE])
storage.mode(count_mat) <- "integer"

# Peak coordinates for BED output
peak_coords <- data.frame(
  peakID = fc$Geneid,
  chr    = fc$Chr,
  start  = as.integer(fc$Start) - 1L,   # convert back to 0-based BED
  end    = as.integer(fc$End),
  stringsAsFactors = FALSE
)
rownames(peak_coords) <- peak_coords$peakID

# --- 2. Load hkTSS scale factors as DESeq2 size factors ---------------------
# hkTSS_scale_factors.tsv columns: sample, hk_tss_reads, hkTSS_scaleFactor
# DESeq2 sizeFactor = 1 / scaleFactor
hk_sf      <- read.table(hk_sf_file, header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE)
hk_sf_vec  <- 1 / hk_sf$hkTSS_scaleFactor
names(hk_sf_vec) <- hk_sf$sample

# Subset to samples present in count matrix
common_samples <- intersect(names(hk_sf_vec), colnames(count_mat))
count_mat      <- count_mat[, common_samples, drop = FALSE]
samples        <- samples[common_samples, , drop = FALSE]
hk_sf_vec      <- hk_sf_vec[common_samples]

message(sprintf("Samples used: %s", paste(common_samples, collapse = ", ")))

# --- 3. Basic filtering ------------------------------------------------------
keep      <- rowSums(count_mat >= 5) >= 2
count_mat <- count_mat[keep, ]
message(sprintf("Peaks after filtering (>= 5 reads in >= 2 samples): %d", nrow(count_mat)))

# --- 4. Build full DESeq2 object with hkTSS size factors --------------------
dds_full <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = samples,
  design    = ~ condition
)
sizeFactors(dds_full) <- hk_sf_vec[colnames(dds_full)]

# --- 5. Run comparisons: each condition vs DMSO -----------------------------
conditions <- unique(samples$condition)
conditions <- conditions[conditions != "DMSO"]

# Thresholds
LFC_THRESH <- 1.0    # |log2FC| >= 1 (2-fold change)
FDR_THRESH <- 0.05

all_gained <- list()
all_lost   <- list()

for (cond in conditions) {

  message(sprintf("\n-- %s vs DMSO --", cond))

  # Subset to this condition + DMSO
  sel      <- samples$condition %in% c("DMSO", cond)
  sub_mat  <- count_mat[, sel, drop = FALSE]
  sub_meta <- samples[sel, , drop = FALSE]
  sub_meta$condition <- factor(sub_meta$condition,
                                levels = c("DMSO", cond))  # DMSO = reference

  dds <- DESeqDataSetFromMatrix(
    countData = sub_mat,
    colData   = sub_meta,
    design    = ~ condition
  )
  sizeFactors(dds) <- hk_sf_vec[colnames(dds)]

  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds,
                 contrast     = c("condition", cond, "DMSO"),
                 alpha        = FDR_THRESH,
                 lfcThreshold = 0)

  res_df <- as.data.frame(res)
  res_df$peakID <- rownames(res_df)
  res_df <- merge(res_df, peak_coords, by = "peakID")

  # Save full results
  write.csv(res_df,
            file = file.path(out_dir, sprintf("%s_vs_DMSO_DESeq2_results.csv", cond)),
            row.names = FALSE)

  # Gained: significantly more open in treatment vs DMSO
  gained <- res_df[!is.na(res_df$padj) &
                     res_df$padj  < FDR_THRESH &
                     res_df$log2FoldChange >= LFC_THRESH, ]

  # Lost: significantly more closed in treatment vs DMSO
  lost   <- res_df[!is.na(res_df$padj) &
                     res_df$padj  < FDR_THRESH &
                     res_df$log2FoldChange <= -LFC_THRESH, ]

  message(sprintf("  Gained: %d peaks | Lost: %d peaks", nrow(gained), nrow(lost)))

  # Write BED files (sorted)
  write_bed <- function(df, path) {
    if (nrow(df) == 0) {
      message(sprintf("  [!] No peaks for %s - writing empty BED", basename(path)))
      writeLines("", path)
      return(invisible(NULL))
    }
    df <- df[order(df$chr, df$start), ]
    bed <- data.frame(chr   = df$chr,
                      start = df$start,
                      end   = df$end,
                      name  = df$peakID,
                      score = round(-log10(df$padj + 1e-300), 2),
                      strand = ".")
    write.table(bed, file = path, sep = "\t",
                quote = FALSE, row.names = FALSE, col.names = FALSE)
  }

  gained_bed <- file.path(out_dir, sprintf("%s_vs_DMSO_gained.bed", cond))
  lost_bed   <- file.path(out_dir, sprintf("%s_vs_DMSO_lost.bed",   cond))
  write_bed(gained, gained_bed)
  write_bed(lost,   lost_bed)

  all_gained[[cond]] <- gained
  all_lost[[cond]]   <- lost
}

# --- 6. Summary table -------------------------------------------------------
summary_df <- data.frame(
  comparison = paste0(conditions, "_vs_DMSO"),
  gained     = sapply(conditions, function(c) nrow(all_gained[[c]])),
  lost       = sapply(conditions, function(c) nrow(all_lost[[c]]))
)
message("\n-- Summary --")
print(summary_df)
write.csv(summary_df,
          file = file.path(out_dir, "diff_atac_summary.csv"),
          row.names = FALSE)

message("\nDone: diff_atac_vs_dmso.R")
