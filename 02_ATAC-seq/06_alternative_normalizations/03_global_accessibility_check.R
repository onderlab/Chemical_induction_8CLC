#!/usr/bin/env Rscript
# =========================================================
# global_accessibility_check.R
# ---------------------------------------------------------
# Compares global chromatin accessibility before and after
# normalization. For each sample:
#   1. Raw total counts (pre-normalization)
#   2. Total signal under CPM, stablePeak and hkTSS normalizations
#   3. Fold change relative to the DMSO mean + bar plots
#
# Usage:
#   Rscript global_accessibility_check.R \
#     <consensus_peak_counts.txt> \
#     <samples.tsv> \
#     <hkTSS_scale_factors.tsv> \
#     <stablePeak_scale_factors.tsv> \
#     <out_dir>
# =========================================================

args          <- commandArgs(trailingOnly = TRUE)
count_file    <- args[1]
sample_file   <- args[2]
hk_sf_file    <- args[3]
sp_sf_file    <- args[4]
out_dir       <- args[5]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(tidyr)
  library(dplyr)
})

# --- 1. Load data -----------------------------------------------------------
samples <- as.data.frame(fread(sample_file))
samples <- samples[samples[[1]] != "sample", ]
colnames(samples)[1:2] <- c("sample", "condition")

fc <- as.data.frame(fread(count_file, skip = 1))
count_cols <- grep("\\.shifted\\.bam$", colnames(fc), value = TRUE)
count_mat  <- fc[, count_cols, drop = FALSE]
rownames(count_mat) <- fc$Geneid
colnames(count_mat) <- sub("\\.shifted\\.bam$", "", basename(count_cols))
count_mat  <- as.matrix(count_mat[, samples$sample, drop = FALSE])
storage.mode(count_mat) <- "integer"

# --- 2. Load scale factors --------------------------------------------------
hk_sf <- read.table(hk_sf_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sp_sf <- read.table(sp_sf_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# --- 3. Compute global signal per sample ------------------------------------
# a) Raw total counts (pre-normalization)
total_raw <- colSums(count_mat)

# b) stablePeak normalized: raw * scaleFactor
sp_scale <- setNames(sp_sf$scaleFactor, sp_sf$sample)
total_sp <- total_raw * sp_scale[names(total_raw)]

# c) hkTSS normalized: raw * hkTSS_scaleFactor
hk_scale <- setNames(hk_sf$hkTSS_scaleFactor, hk_sf$sample)
total_hk <- total_raw * hk_scale[names(total_raw)]

# d) CPM: total / lib_size * 1e6, then sum
lib_sizes <- colSums(count_mat)
cpm_mat   <- t(t(count_mat) / lib_sizes) * 1e6
total_cpm <- colSums(cpm_mat)

# --- 4. Combine -------------------------------------------------------------
df <- data.frame(
  sample     = names(total_raw),
  condition  = samples$condition[match(names(total_raw), samples$sample)],
  raw        = as.numeric(total_raw),
  CPM        = as.numeric(total_cpm),
  stablePeak = as.numeric(total_sp),
  hkTSS      = as.numeric(total_hk),
  stringsAsFactors = FALSE
)

# Normalize to the DMSO mean (fold change)
dmso_mean_raw <- mean(df$raw[df$condition == "DMSO"])
dmso_mean_cpm <- mean(df$CPM[df$condition == "DMSO"])
dmso_mean_sp  <- mean(df$stablePeak[df$condition == "DMSO"])
dmso_mean_hk  <- mean(df$hkTSS[df$condition == "DMSO"])

df$raw_fc        <- df$raw        / dmso_mean_raw
df$CPM_fc        <- df$CPM        / dmso_mean_cpm
df$stablePeak_fc <- df$stablePeak / dmso_mean_sp
df$hkTSS_fc      <- df$hkTSS      / dmso_mean_hk

# --- 5. Write ---------------------------------------------------------------
write.csv(df, file.path(out_dir, "global_accessibility_summary.csv"), row.names = FALSE)

message("Global accessibility (fold change vs DMSO mean):")
print(df[, c("sample","condition","raw_fc","CPM_fc","stablePeak_fc","hkTSS_fc")])

# --- 6. Plot: fold-change comparison ----------------------------------------
fc_long <- df %>%
  select(sample, condition, raw_fc, stablePeak_fc, hkTSS_fc) %>%
  pivot_longer(cols = c(raw_fc, stablePeak_fc, hkTSS_fc),
               names_to = "method", values_to = "fold_change") %>%
  mutate(method = recode(method,
    raw_fc        = "Raw counts",
    stablePeak_fc = "stablePeak",
    hkTSS_fc      = "hkTSS"
  ))

fc_long$method    <- factor(fc_long$method, levels = c("Raw counts","stablePeak","hkTSS"))
fc_long$condition <- factor(fc_long$condition,
                             levels = c("DMSO","PY60","PRBJN","PRBJNsorted","DUX4sorted"))

cond_colors <- c(
  DMSO        = "#888780",
  PY60        = "#1D9E75",
  PRBJN       = "#BA7517",
  PRBJNsorted = "#D4537E",
  DUX4sorted  = "#7F77DD"
)

p1 <- ggplot(fc_long, aes(x = sample, y = fold_change, fill = condition)) +
  geom_col(show.legend = TRUE) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40", linewidth = 0.6) +
  facet_wrap(~method, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = cond_colors) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7)) +
  labs(x = NULL, y = "Fold change vs DMSO",
       title = "Global chromatin accessibility - fold change vs DMSO",
       subtitle = "Values > 1 indicate global increase in accessibility",
       fill = "Condition")

ggsave(file.path(out_dir, "global_accessibility_foldchange.pdf"),
       p1, width = 16, height = 6)

# --- 7. Plot: raw total counts (pre-normalization) --------------------------
p2 <- ggplot(df, aes(x = sample, y = raw / 1e6, fill = condition)) +
  geom_col(show.legend = TRUE) +
  scale_fill_manual(values = cond_colors) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  labs(x = NULL, y = "Total reads in peaks (millions)",
       title = "Total reads in consensus peaks - before normalization",
       fill = "Condition")

ggsave(file.path(out_dir, "global_accessibility_raw_counts.pdf"),
       p2, width = 10, height = 6)

message(sprintf("\nDone. Files: %s", out_dir))
message("  global_accessibility_summary.csv")
message("  global_accessibility_foldchange.pdf")
message("  global_accessibility_raw_counts.pdf")
