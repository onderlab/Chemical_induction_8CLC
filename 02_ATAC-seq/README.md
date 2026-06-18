# ATAC-seq Analysis

End-to-end ATAC-seq workflow from aligned BAMs to normalized signal tracks,
differential accessibility, embryo comparison, motif enrichment and TF
footprinting.

All file paths are placeholders (`path/to/...`) and must be adapted to your
own environment. The pipeline was developed on macOS; the `sed -i ''` calls
and the `ulimit -n` line are macOS-specific (Linux equivalents are included
as fallbacks where relevant).

---

## Folder layout

```
08_ATAC-seq/
├── 01_main_pipeline/
│   ├── 01_run_Atac_from_bam.sh        # 12-step master pipeline
│   ├── 02_atac_multinorm.R            # stable-peak size factors + PCA
│   ├── 03_hkTSS_normalize.R           # housekeeping-TSS scale factors
│   └── 04_compare_region_sets.R       # normalization comparison summaries
│
├── 02_differential_analysis/
│   ├── 01_diff_atac_vs_dmso.R         # DESeq2 each condition vs DMSO -> gained/lost BEDs
│   └── 02_all_conditions_heatmaps.sh  # gained/lost heatmaps (hkTSS bigWigs)
│
├── 03_embryo_comparison/
│   ├── 01_embryo_liftover_gained_lost_heatmaps.sh  # GSE101571 8C/ICM hg19->hg38 + heatmaps
│   ├── 02_zga_embryo_heatmap.sh                     # embryo signal over ZGA TSS
│   └── 03_8c_geneset_heatmap.sh                     # embryo + our samples over 8C geneset
│
├── 04_motif_analysis/
│   ├── 01_motif_enrichment_gained_FIMO_AME.sh  # MEME Suite (AME + FIMO)
│   ├── 02_homer_motif_enrichment.sh            # HOMER known + de novo motifs
│   └── 03_tprx1_footprint_motif.sh             # FIMO scan at a TPRX1 footprint
│
├── 05_footprinting_TOBIAS/
│   └── 01_tobias_parallel.sh          # ATACorrect -> ScoreBigwig -> BINDetect
│
└── 06_alternative_normalizations/
    ├── 01_run_downsampled_pipeline.sh # depth-matched control pipeline
    ├── 02_rip_normalize_atac.sh       # reads-in-peaks normalization
    ├── 03_global_accessibility_check.R# global accessibility summary table
    └── 04_plot_global_accessibility.R # publication plots from that summary
```

---

## Execution order

1. `01_main_pipeline/01_run_Atac_from_bam.sh` runs the full core workflow
   (Steps 1-12) and calls the three R scripts in the same folder. It expects:
   - `samples.tsv` with columns: `sample`, `condition`, `replicate`, `bam`
   - `zga_genes.txt` and `housekeeping_genes.txt` (one gene symbol per line)
   - a reference genome directory with GTF, chrom.sizes and ENCODE blacklist
2. `02_differential_analysis` for DESeq2 gained/lost peaks and their heatmaps.
3. `03_embryo_comparison` to overlay published 8C/ICM embryo ATAC signal.
4. `04_motif_analysis` and `05_footprinting_TOBIAS` for TF analysis.
5. `06_alternative_normalizations` are control/robustness analyses
   (downsampling, reads-in-peaks, global accessibility quantification).

---

## Normalization strategies

This dataset includes conditions (e.g. DUX4-driven) that undergo genome-wide
chromatin opening, which violates the equal-loading assumption of standard
library-size normalization. Several normalizations are therefore provided and
compared:

- **CPM / RPGC** - standard library-size normalizations (deepTools).
- **stablePeak** - DESeq2 size factors estimated only from condition-invariant
  "stable" peaks (`atac_multinorm.R`), robust to global shifts.
- **hkTSS** - scale factors equalising read counts at housekeeping-gene TSS
  windows (`hkTSS_normalize.R`); used as size factors in the differential
  analysis.
- **RIP** - reads-in-peaks normalization (`rip_normalize_atac.sh`).
- **downsampling** - a depth-matched control to rule out library-depth artefacts.

`compare_region_sets.R` and `global_accessibility_check.R` quantify how each
method behaves (e.g. housekeeping-TSS CV across samples).

---

## Software requirements

Command-line: samtools, bedtools, deepTools (alignmentSieve, bamCoverage,
computeMatrix, plotProfile, plotHeatmap), MACS2, featureCounts (Subread),
MEME Suite (ame, fimo, fasta-get-markov), HOMER, TOBIAS, UCSC tools
(liftOver, wigToBigWig, bigWigToBedGraph, bedGraphToBigWig), gawk, Python 3
(pandas).

R: DESeq2, data.table, matrixStats, ggplot2, Rsamtools, GenomicRanges,
GenomicAlignments, tidyr, dplyr.

---

## Published data used

- **GSE101571** - human 8-cell and ICM ATAC-seq (hg19; lifted over to hg38).
- **JASPAR2024** CORE vertebrates - motif database for FIMO / AME / TOBIAS.
- ENCODE hg38 blacklist v2.
