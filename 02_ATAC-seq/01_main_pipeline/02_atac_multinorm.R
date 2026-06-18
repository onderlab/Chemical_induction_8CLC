#!/usr/bin/env Rscript
# =========================================================
# atac_multinorm.R
# ---------------------------------------------------------
# Stable-peak based size-factor estimation for ATAC-seq.
#   - Selects "stable" peaks: high mean CPM and low CV across samples
#     (interpreted as condition-invariant technical references)
#   - Estimates DESeq2 size factors using these as controlGenes
#   - Provides QC on housekeeping-promoter peaks
#   - Outputs a stablePeak_scale_factors.tsv to be passed to
#     bamCoverage --scaleFactor in downstream steps
# =========================================================

args <- commandArgs(trailingOnly = TRUE)
count_file  <- args[1]
sample_file <- args[2]
hk_bed      <- args[3]
out_dir     <- args[4]
qc_dir      <- args[5]

suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(matrixStats)
  library(ggplot2)
})

# --- Read sample info ---
samples <- as.data.frame(fread(sample_file))
rownames(samples) <- samples$sample

# --- Read featureCounts output ---
fc <- as.data.frame(fread(count_file, skip = 1))
count_cols <- grep("\\.shifted\\.bam$", colnames(fc), value = TRUE)
count_mat <- fc[, count_cols, drop = FALSE]
rownames(count_mat) <- fc$Geneid
colnames(count_mat) <- sub("\\.shifted\\.bam$", "", basename(count_cols))
count_mat <- as.matrix(count_mat[, samples$sample, drop = FALSE])
storage.mode(count_mat) <- "integer"

# --- Peak coordinates for HK overlap ---
peak_coords <- data.frame(
  peakID = fc$Geneid,
  chr    = fc$Chr,
  start  = fc$Start,
  end    = fc$End,
  stringsAsFactors = FALSE
)

# --- Read HK TSS BED and find overlapping peaks ---
hk_tss <- as.data.frame(fread(hk_bed, header = FALSE))
colnames(hk_tss) <- c("chr", "start", "end", "gene", "score", "strand")
hk_tss$prom_start <- pmax(0, hk_tss$start - 2000)
hk_tss$prom_end   <- hk_tss$end + 2000

hk_peak_idx <- sapply(seq_len(nrow(peak_coords)), function(i) {
  any(peak_coords$chr[i] == hk_tss$chr &
      peak_coords$start[i] < hk_tss$prom_end &
      peak_coords$end[i]   > hk_tss$prom_start)
})
hk_peaks <- peak_coords$peakID[hk_peak_idx]
message(sprintf("HK promoter-overlapping peaks: %d / %d total",
                length(hk_peaks), nrow(peak_coords)))

# --- Initial filtering ---
dds0 <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = samples,
  design    = ~ condition
)
keep <- rowSums(counts(dds0) >= 10) >= 2
dds0 <- dds0[keep, ]
raw_counts <- counts(dds0)

# --- Stable-peak selection (CPM-based mean and CV filter) ---
# Stable peaks: not low-expressed (top 60% by mean CPM) AND
#               low variability (bottom 20% by CV across samples)
lib_sizes <- colSums(raw_counts)
cpm <- t(t(raw_counts) / lib_sizes) * 1e6
peak_mean <- rowMeans(cpm)
peak_cv   <- rowSds(as.matrix(cpm)) / (peak_mean + 1e-8)

stable_idx <- which(
  peak_mean >= quantile(peak_mean, 0.40) &
  peak_cv   <= quantile(peak_cv,   0.20)
)
stable_peaks <- rownames(raw_counts)[stable_idx]

message(sprintf("Stable peaks selected: %d", length(stable_peaks)))

hk_in_stable <- intersect(stable_peaks, hk_peaks)
message(sprintf("HK promoter peaks among stable peaks: %d / %d HK peaks",
                length(hk_in_stable), length(hk_peaks)))

write.table(stable_peaks, file = file.path(out_dir, "stable_peaks.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# --- DESeq2 size factors estimated from stable peaks only ---
# estimateSizeFactors with controlGenes uses only those rows as reference,
# avoiding the influence of condition-driven peak changes on normalization.
dds <- DESeqDataSetFromMatrix(
  countData = raw_counts,
  colData   = samples,
  design    = ~ condition
)
dds <- estimateSizeFactors(dds, controlGenes = which(rownames(dds) %in% stable_peaks))
sf <- sizeFactors(dds)
scale_factor <- 1 / sf

sf_tab <- data.frame(
  sample      = names(sf),
  sizeFactor  = as.numeric(sf),
  scaleFactor = as.numeric(scale_factor)
)
write.table(sf_tab, file = file.path(out_dir, "stablePeak_scale_factors.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- HK promoter normalization QC ---
hk_in_counts <- intersect(hk_peaks, rownames(raw_counts))
if (length(hk_in_counts) >= 10) {
  hk_cpm <- cpm[hk_in_counts, , drop = FALSE]
  hk_mean_per_sample <- colMeans(hk_cpm)
  hk_qc <- data.frame(
    sample           = names(hk_mean_per_sample),
    hk_mean_cpm      = as.numeric(hk_mean_per_sample),
    stablePeak_sf    = as.numeric(sf[names(hk_mean_per_sample)]),
    hk_corrected_cpm = as.numeric(hk_mean_per_sample *
                                  scale_factor[names(hk_mean_per_sample)])
  )
  write.csv(hk_qc, file = file.path(qc_dir, "HK_promoter_normalization_QC.csv"),
            row.names = FALSE)
  message("HK promoter QC table written.")
} else {
  message("Too few HK peaks in count matrix for QC table.")
}

# --- PCA on VST-normalized data ---
vsd <- vst(dds, blind = TRUE)
pca <- prcomp(t(assay(vsd)))
pca_df <- data.frame(
  sample    = rownames(pca$x),
  PC1       = pca$x[, 1],
  PC2       = pca$x[, 2],
  condition = samples[rownames(pca$x), "condition"]
)
percentVar <- (pca$sdev^2) / sum(pca$sdev^2)

p <- ggplot(pca_df, aes(PC1, PC2, color = condition, label = sample)) +
  geom_point(size = 4) +
  geom_text(vjust = -0.7, size = 3) +
  theme_bw(base_size = 12) +
  xlab(sprintf("PC1: %.1f%%", percentVar[1] * 100)) +
  ylab(sprintf("PC2: %.1f%%", percentVar[2] * 100))

ggsave(file.path(qc_dir, "PCA_stablePeakNorm.pdf"), p, width = 8, height = 6)
message("Done: atac_multinorm.R")
