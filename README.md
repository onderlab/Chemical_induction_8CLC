# RNA-seq, TE, and scRNA-seq Analysis Code

This repository contains the analysis scripts used in the manuscript.
All file paths in the scripts are placeholders (e.g. `path/to/...`) and must be
adapted to your local environment before running.

---

## Repository structure

```
.
├── 01_RNA-seq/
│   ├── 01_salmon_quantification_and_DESeq2.md   # Salmon + DESeq2 pipeline
│   └── 02_quant_to_gene_TPM_matrix.R            # Transcript -> gene-level TPM
│
├── 02_TE_analysis/
│   ├── 01_SalmonTE_pipeline.md                  # SalmonTE-based TE DE analysis
│   ├── 02_TEtranscripts_pipeline.md             # STAR + TEtranscripts pipeline
│   └── 03_TE_heatmap_generation.R               # TE heatmap from TPM values
│
├── 03_scRNA-seq/
│   └── 01_embryo_atlas_integration_and_mapping.R
│       # Seurat integration of Yan/Petropoulos/Mazid + MapQuery of our data
│
├── 04_GSEA/
│   ├── 01_GSEA_fgsea_MSigDB_collections.R       # fgsea on Hallmark/Reactome/GO BP
│   ├── 02_GSEA_clusterProfiler_custom_gene_sets.R
│   ├── 03_GSEA_NES_vs_FDR_plot.R                # NES vs FDR scatter plot
│   └── 04_GSEA_custom_geneset_from_DESeq2.R     # GSEA from DESeq2 stat ranking
│
├── 05_Heatmaps/
│   ├── 01_heatmap_VST_normalization_tximport.R  # VST heatmap (tximport input)
│   ├── 02_embryo_heatmap_TPM_zscore.R           # Embryo-stage averaged heatmap
│   └── 03_ordered_geneset_heatmap.R             # Gene-set-ordered heatmap
│
├── 06_Volcano_plots/
│   └── 01_volcano_from_DESeq2_results.R         # Volcano from DESeq2 results
│
├── 07_Correlation_analysis/
│   └── 01_DUX4_targets_pearson_correlation.R    # DUX4 target gene correlation
│
└── 08_ATAC-seq/                                 # full ATAC-seq workflow (see its own README)
    ├── 01_main_pipeline/                         # BAM -> peaks -> normalized tracks (12 steps)
    ├── 02_differential_analysis/                 # DESeq2 gained/lost peaks + heatmaps
    ├── 03_embryo_comparison/                     # 8C/ICM embryo ATAC overlay
    ├── 04_motif_analysis/                        # MEME Suite + HOMER motif enrichment
    ├── 05_footprinting_TOBIAS/                   # TF footprinting (ATACorrect/BINDetect)
    └── 06_alternative_normalizations/            # downsampling / RIP / global accessibility
```

---


## Software / package requirements

### Command-line
- FastQC
- Salmon (v1.1.0 or later)
- STAR
- samtools
- TEtranscripts
- SalmonTE
- Python 3 (pandas)
- **ATAC-seq:** bedtools, deepTools, MACS2, featureCounts (Subread),
  MEME Suite, HOMER, TOBIAS, UCSC tools (liftOver, wigToBigWig, etc.), gawk

### R packages
- **Core:** `tximport`, `DESeq2`, `GenomicFeatures`, `biomaRt`
- **Visualization:** `ggplot2`, `EnhancedVolcano`, `pheatmap`, `ComplexHeatmap`,
  `circlize`, `RColorBrewer`, `ggrepel`
- **Enrichment:** `fgsea`, `msigdbr`, `clusterProfiler`, `enrichplot`, `GSEABase`
- **scRNA-seq:** `Seurat`, `dplyr`
- **I/O:** `readxl`, `openxlsx`, `writexl`, `readr`, `data.table`

---

## Reference datasets used in the scRNA-seq atlas

- **Yan et al. (2013)** -- *Nat. Struct. Mol. Biol.*, PMID: 23934149
- **Petropoulos et al. (2016)** -- *Cell*, PMID: 27062923
- **Mazid et al. (2022)** -- *Nature*, PMID: 35314832

Integration was performed using Seurat's anchor-based integration
(Stuart, Butler et al., *Cell*, 2019).

---

## Notes

- All `path/to/...` placeholders in the scripts must be replaced with the
  appropriate paths on your system.
- Several scripts assume an Excel file as input for DESeq2 results or gene sets;
  these were exported from earlier steps in the pipeline.

