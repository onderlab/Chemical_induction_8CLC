# Salmon Quantification and DESeq2 Differential Expression Analysis

## 1. Quality Control with FastQC

```bash
fastqc --threads 6 *.fastq.gz
```

## 2. Download Salmon Index for Target Transcriptome

Salmon index assets are available at: http://refgenomes.databio.org/

Download the hg38 transcriptome index for Salmon (produced with `salmon index` using the
selective alignment method, which improves quantification accuracy compared to the
regular index):

`salmon_sa_index:default`

```
http://refgenomes.databio.org/v2/asset/hg38/salmon_sa_index/archive?tag=default
```

## 3. Run Salmon for Each Sample

Use `--validateMappings`, `--seqBias`, and `--gcBias` parameters. Add the location of the
index folder, etc.:

```bash
salmon quant \
  --threads 8 \
  --validateMappings \
  --seqBias \
  --gcBias \
  -i path/to/salmon_sa_index \
  --libType A \
  -1 path/to/fastq/Sample_1.fastq.gz \
  -2 path/to/fastq/Sample_2.fastq.gz \
  -o path/to/output/quant/Sample
```

## 4. Download Genome/Transcriptome Annotation File

Download the GENCODE primary assembly annotation GTF from:
https://www.gencodegenes.org/human/

For example: `gencode.v35.primary_assembly.annotation.gtf`

```
ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.primary_assembly.annotation.gtf.gz
```

## 5. DESeq2 Analysis in R

Put Salmon's output folder named `quant` in your working directory, and include the
`gencode.v35.primary_assembly.annotation.gtf` file.

```r
# ---- Install required packages ----
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("tximport")
BiocManager::install("GenomicFeatures")
install.packages(c("Rcpp", "readr"))

# ---- Load libraries ----
library(tximport)
library(GenomicFeatures)
library(readr)
library(DESeq2)
library(gplots)
library(ggplot2)
library(RColorBrewer)
library(psych)
library(EnhancedVolcano)

setwd("path/to/working/directory")

# ---- Build tx2gene mapping from GTF ----
txdb <- makeTxDbFromGFF("gencode.v35.primary_assembly.annotation.gtf")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, keys = k, keytype = "TXNAME", columns = "GENEID")

head(tx2gene)

# ---- Import Salmon quantifications ----
# NOTE: countsFromAbundance = "no" (the default) passes raw counts plus an
# averageTranscriptLength offset to DESeq2. This is the correct setting for
# differential expression analysis - DESeq2 itself models the count
# distribution and applies length correction internally.
samples <- read.table("samples.txt", header = TRUE)
files   <- file.path("quant", samples$sample, "quant.sf")
names(files) <- paste0(samples$sample)

txi.salmon <- tximport(files,
                       type    = "salmon",
                       tx2gene = tx2gene,
                       countsFromAbundance = "no")

# ---- Inspect ----
head(txi.salmon$counts)

# ---- Export raw count matrix ----
write.table(txi.salmon$counts, "txi_salmon_counts_ENSG.txt", sep = "\t")
write.table((txi.salmon$counts[, 0]), "ENSG_IDs.txt", sep = "\t")

# ---- Prepare gene names from Ensembl BioMart in order of ENSG IDs ----
genenames <- as.matrix(read.table("GeneNames.txt", header = FALSE))

# ============================================================
# DESeq2 Analysis
# ============================================================
# Set the reference level explicitly so fold changes have a known direction.
samples$condition <- relevel(factor(samples$condition), ref = "DMSO")

dds <- DESeqDataSetFromTximport(txi.salmon, samples, ~ condition)
rownames(dds) <- genenames

# ---- Group-size-aware prefiltering ----
# Keep genes with >=10 counts in at least as many samples as the
# smallest group (recommended in the DESeq2 vignette).
smallestGroupSize <- min(table(samples$condition))
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds  <- dds[keep, ]

# ---- Run DESeq2 ----
dds <- DESeq(dds)

# ---- Use explicit contrasts ----
# results() with no contrast uses the last factor level vs. the first,
# which can be ambiguous. Always specify direction:
#   contrast = c("condition", numerator, denominator)
res <- results(dds, contrast = c("condition", "PY60_Sorted", "DMSO"))

write.table(res, "DESeq2Results_PY60_Sorted_vs_DMSO.txt", sep = "\t")
write.table(counts(dds), "DESeq2CountsFiltered.txt", sep = "\t")

summary(res)
```
