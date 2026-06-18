#!/usr/bin/env Rscript
# =========================================================
# compare_region_sets.R
# ---------------------------------------------------------

# =========================================================

args <- commandArgs(trailingOnly = TRUE)
matrix_dir  <- args[1]
sample_file <- args[2]
out_dir     <- args[3]

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

samples <- as.data.frame(fread(sample_file))
modes   <- c("CPM", "RPGC", "stablePeak")
regions <- c("allTSS", "hkTSS", "zgaTSS")

all_res <- list()

for (mode in modes) {
  for (reg in regions) {
    tsv_file <- file.path(matrix_dir, paste0(reg, "_", mode, "_profile.tsv"))
    if (!file.exists(tsv_file)) {
      message(sprintf("Missing: %s - skipping", tsv_file))
      next
    }

    # plotProfile --outFileNameData TSV format:
    #   Line 1: "bin labels\t<label1>\t<label2>..." (one label per sample group)
    #   Line 2+: "<sample>\tgenes\t<val1>\t<val2>..."  (one row per sample)
    # fread() cannot handle this because the header and data rows have
    # different column counts, so parse manually.
    lines <- readLines(tsv_file)

    header_idx <- grep("^bin labels", lines, ignore.case = TRUE)
    if (length(header_idx) == 0) {
      message(sprintf("Could not find 'bin labels' header in: %s - skipping", tsv_file))
      next
    }

    data_lines <- lines[(header_idx[1] + 1):length(lines)]
    data_lines <- data_lines[nchar(trimws(data_lines)) > 0]  # drop blank lines

    # Each data line: sample_name \t "genes" \t val1 \t val2 ...
    parsed <- lapply(data_lines, function(ln) {
      parts <- strsplit(ln, "\t")[[1]]
      sample_name <- parts[1]
      values      <- suppressWarnings(as.numeric(parts[-(1:2)]))  # skip sample + "genes"
      list(sample = sample_name, values = values)
    })

    n_bins <- length(parsed[[1]]$values)
    if (n_bins == 0) {
      message(sprintf("No numeric bin values found in: %s - skipping", tsv_file))
      next
    }

    num_mat      <- do.call(rbind, lapply(parsed, `[[`, "values"))
    sample_names <- sapply(parsed, `[[`, "sample")

    # Mean signal across all bins per sample
    mean_signal <- rowMeans(num_mat, na.rm = TRUE)
    # Signal at TSS (center bin)
    center_idx  <- ceiling(ncol(num_mat) / 2)
    tss_signal  <- num_mat[, center_idx]

    res <- data.frame(
      sample      = sample_names,
      region      = reg,
      norm_method = mode,
      mean_signal = mean_signal,
      tss_signal  = tss_signal,
      stringsAsFactors = FALSE
    )
    all_res[[paste(mode, reg, sep = "_")]] <- res
  }
}

if (length(all_res) == 0) {
  message("No profile TSV files found. Run the pipeline first.")
  quit(status = 0)
}

summary_df <- do.call(rbind, all_res)
rownames(summary_df) <- NULL

write.csv(summary_df,
          file = file.path(out_dir, "normalization_comparison_summary.csv"),
          row.names = FALSE)

# --- Barplot: TSS signal per sample, faceted by region and norm method ---
p <- ggplot(summary_df, aes(x = sample, y = tss_signal, fill = sample)) +
  geom_col(show.legend = FALSE) +
  facet_grid(region ~ norm_method, scales = "free_y") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7)) +
  labs(x = NULL, y = "Signal at TSS",
       title = "TSS signal comparison across normalizations")

ggsave(file.path(out_dir, "normalization_TSS_signal_comparison.pdf"),
       p, width = 14, height = 10)

# --- HK should be flat across samples: compute CV per norm method ---
hk_rows <- summary_df[summary_df$region == "hkTSS", ]
hk_cv <- aggregate(tss_signal ~ norm_method, data = hk_rows,
                   FUN = function(x) sd(x) / mean(x))
colnames(hk_cv)[2] <- "CV_across_samples"

message("\nHousekeeping TSS signal CV across samples (lower = better normalization):")
print(hk_cv)

write.csv(hk_cv,
          file = file.path(out_dir, "HK_TSS_CV_by_normalization.csv"),
          row.names = FALSE)

# --- ZGA enrichment ratio: ZGA signal / HK signal per sample ---
zga_rows <- summary_df[summary_df$region == "zgaTSS", ]
merged <- merge(
  zga_rows[, c("sample", "norm_method", "tss_signal")],
  hk_rows[, c("sample", "norm_method", "tss_signal")],
  by = c("sample", "norm_method"),
  suffixes = c("_zga", "_hk")
)
merged$zga_over_hk <- merged$tss_signal_zga / (merged$tss_signal_hk + 1e-8)

write.csv(merged,
          file = file.path(out_dir, "ZGA_vs_HK_enrichment_ratio.csv"),
          row.names = FALSE)

p2 <- ggplot(merged, aes(x = sample, y = zga_over_hk, fill = sample)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~ norm_method, scales = "free_y") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7)) +
  labs(x = NULL, y = "ZGA / HK signal ratio",
       title = "ZGA enrichment relative to housekeeping (per normalization)")

ggsave(file.path(out_dir, "ZGA_vs_HK_enrichment_ratio.pdf"),
       p2, width = 12, height = 6)

message("Done: compare_region_sets.R")
