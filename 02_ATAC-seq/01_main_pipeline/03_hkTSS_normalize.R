#!/usr/bin/env Rscript
# =========================================================
# hkTSS_normalize.R
# ---------------------------------------------------------
# Reference-region (HK TSS) based scale-factor computation.
#
# Logic:
#   1. For every sample, count raw paired-end reads falling in
#      housekeeping gene TSS +/- 500 bp windows (from shifted BAM).
#   2. The sample with the LOWEST read count is taken as the reference.
#   3. scale_factor = min_reads / sample_reads
#      -> reference gets 1.0, others get < 1.0 (scaled DOWN to match).
#   4. Pass scale_factor to bamCoverage --scaleFactor
#      (with --normalizeUsing None) to produce hkTSS-normalized bigWigs.
#
# Rationale: housekeeping promoters should be constitutively open and
# unaffected by treatment. Equalising their read counts across samples
# corrects for global accessibility differences that pure library-size
# normalisation cannot.
#
# Usage:
#   Rscript hkTSS_normalize.R \
#     <hk_tss_1bp.bed>      \
#     <shifted_bam_dir>     \
#     <samples.tsv>         \
#     <out_dir>             \
#     [qc_dir]
# =========================================================

args         <- commandArgs(trailingOnly = TRUE)
hk_bed       <- args[1]   # hk_tss_1bp.bed
bam_dir      <- args[2]   # directory containing shifted BAMs
sample_file  <- args[3]   # samples.tsv
out_dir      <- args[4]
qc_dir       <- ifelse(length(args) >= 5, args[5], args[4])

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(qc_dir,  showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(Rsamtools)
  library(GenomicRanges)
  library(GenomicAlignments)
  library(ggplot2)
  library(tidyr)
})

# --- 1. Read HK TSS regions and expand by +/- 500 bp -----------------------
message("Reading HK TSS regions...")
hk <- as.data.frame(fread(hk_bed, header = FALSE))
colnames(hk)[1:4] <- c("chr", "start", "end", "gene")

# +/- 500 bp window (0-based -> GRanges expects 1-based: start+1)
hk_gr <- GRanges(
  seqnames = hk$chr,
  ranges   = IRanges(
    start = pmax(1, hk$start - 500 + 1),
    end   = hk$end + 500
  )
)
hk_gr <- reduce(hk_gr)   # merge overlapping windows
message(sprintf("  HK TSS regions: %d (merged)", length(hk_gr)))

# --- 2. Count reads falling in HK regions per sample -----------------------
samples <- as.data.frame(fread(sample_file))
# skip header row if present
samples <- samples[samples[[1]] != "sample", ]

message("\nHK TSS +/- 500 bp read counts:")
read_counts <- setNames(numeric(nrow(samples)), samples[[1]])

for (i in seq_len(nrow(samples))) {
  sname <- samples[[1]][i]
  bam   <- file.path(bam_dir, paste0(sname, ".shifted.bam"))

  if (!file.exists(bam)) {
    stop(sprintf("BAM not found: %s", bam))
  }

  # Count properly paired reads overlapping HK windows
  param <- ScanBamParam(which = hk_gr, flag = scanBamFlag(isProperPair = TRUE))
  cnt   <- countBam(bam, param = param)
  total <- sum(cnt$records)

  read_counts[sname] <- total
  message(sprintf("  %-25s %d reads", sname, total))
}

# --- 3. Compute scale factors ----------------------------------------------
min_reads  <- min(read_counts)
ref_sample <- names(read_counts)[which.min(read_counts)]
scale_fac  <- min_reads / read_counts

message(sprintf("\nReference sample (lowest reads): %s (%d reads)",
                ref_sample, min_reads))
message("\nScale factors:")
for (s in names(scale_fac)) {
  message(sprintf("  %-25s %.6f  (%d -> ~%d after scaling)",
                  s, scale_fac[s],
                  read_counts[s], round(read_counts[s] * scale_fac[s])))
}

# --- 4. Save scale factors --------------------------------------------------
out_df <- data.frame(
  sample            = names(scale_fac),
  hk_tss_reads      = as.integer(read_counts),
  hkTSS_scaleFactor = as.numeric(scale_fac),
  stringsAsFactors  = FALSE
)

out_tsv <- file.path(out_dir, "hkTSS_scale_factors.tsv")
write.table(out_df, file = out_tsv, sep = "\t",
            quote = FALSE, row.names = FALSE)
message(sprintf("\nScale factors written: %s", out_tsv))

# --- 5. QC plot -------------------------------------------------------------
qc_df <- data.frame(
  sample        = names(read_counts),
  raw_reads     = as.numeric(read_counts),
  scaled_reads  = as.numeric(read_counts * scale_fac),
  stringsAsFactors = FALSE
)
qc_df$sample <- factor(qc_df$sample, levels = qc_df$sample)

qc_long       <- pivot_longer(qc_df, cols = c("raw_reads", "scaled_reads"),
                               names_to = "stage", values_to = "reads")
qc_long$stage <- factor(qc_long$stage,
                         levels = c("raw_reads", "scaled_reads"),
                         labels = c("Raw read counts", "After scaling (equalised)"))

p <- ggplot(qc_long, aes(x = sample, y = reads / 1e6, fill = sample)) +
  geom_col(show.legend = FALSE) +
  geom_hline(yintercept = min_reads / 1e6, linetype = "dashed",
             colour = "grey30", linewidth = 0.6) +
  facet_wrap(~stage, scales = "free_y") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  labs(x = NULL, y = "HK TSS reads (millions)",
       title = "HK TSS +/- 500 bp normalization",
       subtitle = sprintf("Reference: %s (%s million reads)",
                          ref_sample,
                          format(round(min_reads / 1e6, 2), nsmall = 2)))

out_pdf <- file.path(qc_dir, "hkTSS_normalization_QC.pdf")
ggsave(out_pdf, p, width = 14, height = 6)
message(sprintf("QC plot written: %s", out_pdf))

message("\nDone: hkTSS_normalize.R")
