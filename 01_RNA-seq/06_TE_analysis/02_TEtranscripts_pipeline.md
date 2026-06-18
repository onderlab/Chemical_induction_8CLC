# TE Analysis with STAR + TEtranscripts

## 1. STAR alignment (multi-mapping parameters tuned for TE quantification)

```bash
STAR \
  --runThreadN 8 \
  --genomeDir path/to/star/index \
  --readFilesIn path/to/Trimmed_Fastq_Files/Sample_1.trimmed.fq \
                path/to/Trimmed_Fastq_Files/Sample_2.trimmed.fq \
  --outFileNamePrefix path/to/output/Sample/ \
  --winAnchorMultimapNmax 100 \
  --outFilterMultimapNmax  100 \
  --outFilterMismatchNmax  33 \
  --outFilterMismatchNoverLmax 0.3 \
  --alignIntronMin     20 \
  --alignIntronMax     1000000 \
  --alignMatesGapMax   1000000 \
  --alignSJoverhangMin 8 \
  --alignSJDBoverhangMin 1 \
  --outSAMattributes NH HI AS nM XS
```

## 2. Convert SAM to sorted BAM and index

```bash
samtools view -@ 8 -bS path/to/output/Sample/Aligned.out.sam | \
  samtools sort -@ 8 -o path/to/output/Sample/Aligned.sortedByCoord.out.bam

samtools index -@ 8 path/to/output/Sample/Aligned.sortedByCoord.out.bam
```

## 3. Run TEtranscripts for differential TE expression

`--minread 10` removes low-count TEs that would otherwise inflate noise.
`--foldchange 2` corresponds to |log2FC| >= 1, which is a more standard
cut-off than the 1.5 used previously.

```bash
TEtranscripts \
  --mode multi \
  --project TE_DE_PY60_vs_DMSO \
  --stranded no \
  --GTF path/to/gencode.v48.primary_assembly.annotation.gtf \
  --TE  path/to/GRCh38_GENCODE_rmsk_TE.gtf \
  -t path/to/bam/PY60_Sorted_1/Aligned.sortedByCoord.out.bam \
     path/to/bam/PY60_Sorted_2/Aligned.sortedByCoord.out.bam \
  -c path/to/bam/DMSO_1/Aligned.sortedByCoord.out.bam \
     path/to/bam/DMSO_2/Aligned.sortedByCoord.out.bam \
  --DESeq \
  --minread 10 \
  --padj 0.05 \
  --foldchange 2 \
  --sortByFC
```

## 4. Volcano plot from TEtranscripts output (R)

```r
# install.packages(c("data.table","dplyr","ggplot2"))
library(data.table)
library(dplyr)
library(ggplot2)

res_file <- "TE_DE_PY60_vs_DMSO_TEtranscripts_TE_analysis.txt"
res <- fread(res_file)

# TEtranscripts column names are typically:
# id, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj
# (Slight variations possible; verify and adapt if needed.)

# Drop NA padj
res <- res %>% filter(!is.na(padj))

# Thresholds
alpha <- 0.05
lfc   <- 1

sig <- res %>%
  mutate(sig = case_when(
    padj < alpha & log2FoldChange >=  lfc ~ "Up (PY60)",
    padj < alpha & log2FoldChange <= -lfc ~ "Down (PY60)",
    TRUE                                  ~ "NS"
  ))

# Summary
table(sig$sig)

# Volcano plot
ggplot(sig, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(shape = sig), size = 1.4, alpha = 0.7) +
  geom_vline(xintercept = c(-lfc, lfc), linetype = "dashed") +
  geom_hline(yintercept = -log10(alpha),       linetype = "dashed") +
  labs(
    title = "Differentially Expressed Transposable Elements (PY60 vs DMSO)",
    x     = "log2 Fold Change (PY60 / DMSO)",
    y     = "-log10 adjusted p-value"
  ) +
  theme_bw(base_size = 12)
```
