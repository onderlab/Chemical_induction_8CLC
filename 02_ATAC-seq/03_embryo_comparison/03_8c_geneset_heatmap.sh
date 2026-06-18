#!/usr/bin/env bash
# =========================================================
# 8c_geneset_heatmap.sh
# ---------------------------------------------------------
# For the 8C gene set TSS +/- regions, plots:
#   - 8C and ICM embryo ATAC signal
#   - our hkTSS-normalized samples
#
# Prerequisites:
#   - 01_embryo_liftover_gained_lost_heatmaps.sh (hg38 embryo bigWigs)
#   - main pipeline Steps 10-11 (hkTSS bigWigs)
# =========================================================

set -euo pipefail

BASE_DIR="path/to/ATAC_analysis"
PROC_DIR="${BASE_DIR}/processed_multiNorm"
BED_DIR="${PROC_DIR}/beds"
COUNT_DIR="${PROC_DIR}/counts"
BW_DIR="${PROC_DIR}/bigwig"
MATRIX_DIR="${PROC_DIR}/matrices"
HEATMAP_DIR="${PROC_DIR}/heatmaps_diff"
EMBRYO_DIR="${BASE_DIR}/embryo_data"
GENOME_DIR="path/to/hg38"
GTF="${GENOME_DIR}/hg38_genes.gtf"
SAMPLES="${BASE_DIR}/samples.tsv"
THREADS=6

mkdir -p "${BED_DIR}" "${HEATMAP_DIR}"

GENESET_TXT="${BASE_DIR}/8c_geneset.txt"   # gene list (one symbol per line)
BED_8CSET="${BED_DIR}/8c_geneset_tss_1bp.bed"
BW_8C="${EMBRYO_DIR}/GSE101571_8cell_rpkm_hg38.bw"
BW_ICM="${EMBRYO_DIR}/GSE101571_icm_rpkm_hg38.bw"

# --- Step A: 8C gene set -> TSS BED -----------------------------------------
echo "=== A. Building 8C gene set TSS BED ==="

if [[ -s "${BED_8CSET}" ]]; then
    echo "  BED exists - skipping."
else
    echo "  Filtering 8C gene set from all_gene_tss_1bp.bed..."

    ALL_TSS="${BED_DIR}/all_gene_tss_1bp.bed"

    # Build all_gene_tss_1bp.bed from GTF if missing
    if [[ ! -s "${ALL_TSS}" ]]; then
        echo "  Extracting TSS coordinates from GTF..."
        gawk 'BEGIN{FS=OFS="\t"}
        $3=="gene" {
            gene_name=""; strand=$7;
            if (match($9, /gene_name "([^"]+)"/, a)) gene_name=a[1];
            if (gene_name!="") {
                if (strand=="+") {tss=$4-1} else {tss=$5-1}
                print $1,tss,tss+1,gene_name,".",strand
            }
        }' "${GTF}" \
        | awk '$1 ~ /^chr([0-9]+|X|Y)$/' \
        > "${ALL_TSS}"
    fi

    awk 'NR==FNR {keep[$1]=1; next} ($4 in keep)' \
      "${GENESET_TXT}" \
      "${ALL_TSS}" \
      > "${BED_8CSET}"

    N=$(wc -l < "${BED_8CSET}")
    echo "  BED written: ${N} genes matched from gene set"
fi

N_GENES=$(wc -l < "${BED_8CSET}")

# --- Step B: hkTSS bigWig list ----------------------------------------------
BWFILES_HK=()
LABELS_HK=()
while IFS=$'\t' read -r SAMPLE _SF _SC; do
    BWFILES_HK+=("${BW_DIR}/${SAMPLE}.hkTSS.bw")
    LABELS_HK+=("${SAMPLE}")
done < <(tail -n +2 "${COUNT_DIR}/hkTSS_scale_factors.tsv")

# --- Step C: embryo heatmap (8C + ICM) --------------------------------------
echo "=== B. Embryo heatmap (8C + ICM) ==="

MAT_EMB="${MATRIX_DIR}/8c_geneset_embryo.gz"
PDF_EMB="${HEATMAP_DIR}/8c_geneset_embryo_heatmap.pdf"

if [[ ! -f "${MAT_EMB}" ]]; then
    echo "  Computing matrix: embryo"
    computeMatrix reference-point \
      --referencePoint TSS \
      -b 3000 -a 3000 \
      -R "${BED_8CSET}" \
      -S "${BW_8C}" "${BW_ICM}" \
      --skipZeros \
      --missingDataAsZero \
      -p ${THREADS} \
      -o "${MAT_EMB}"
fi

if [[ ! -f "${PDF_EMB}" ]]; then
    echo "  Plotting: embryo heatmap"
    plotHeatmap -m "${MAT_EMB}" \
      --colorMap Reds \
      --whatToShow 'heatmap and colorbar' \
      --samplesLabel "8C" "ICM" \
      --refPointLabel "TSS" \
      --xAxisLabel "Distance from TSS (bp)" \
      --zMin 0 \
      --sortUsing mean \
      --sortRegions descend \
      --missingDataColor 'white' \
      --heatmapHeight 18 \
      --heatmapWidth 3 \
      --plotTitle "8C geneset TSS (n=${N_GENES})" \
      -out "${PDF_EMB}"
    echo "  Saved: ${PDF_EMB}"
fi

# --- Step D: hkTSS-normalized samples heatmap -------------------------------
echo "=== C. hkTSS-normalized samples heatmap ==="

MAT_HK="${MATRIX_DIR}/8c_geneset_hkTSS.gz"
PDF_HK="${HEATMAP_DIR}/8c_geneset_hkTSS_heatmap.pdf"

if [[ ! -f "${MAT_HK}" ]]; then
    echo "  Computing matrix: hkTSS samples"
    computeMatrix reference-point \
      --referencePoint TSS \
      -b 3000 -a 3000 \
      -R "${BED_8CSET}" \
      -S "${BWFILES_HK[@]}" \
      --skipZeros \
      --missingDataAsZero \
      -p ${THREADS} \
      -o "${MAT_HK}"
fi

if [[ ! -f "${PDF_HK}" ]]; then
    echo "  Plotting: hkTSS heatmap"
    plotHeatmap -m "${MAT_HK}" \
      --colorMap Reds \
      --whatToShow 'heatmap and colorbar' \
      --samplesLabel ${LABELS_HK[*]} \
      --refPointLabel "TSS" \
      --xAxisLabel "Distance from TSS (bp)" \
      --zMin 0 \
      --sortUsing mean \
      --sortRegions descend \
      --missingDataColor 'white' \
      --heatmapHeight 18 \
      --heatmapWidth 3 \
      --plotTitle "8C geneset TSS (n=${N_GENES})" \
      -out "${PDF_HK}"
    echo "  Saved: ${PDF_HK}"
fi

echo ""
echo "Done."
echo "  TSS BED       : ${BED_8CSET}"
echo "  Embryo        : ${PDF_EMB}"
echo "  hkTSS samples : ${PDF_HK}"
