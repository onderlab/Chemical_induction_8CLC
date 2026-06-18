#!/usr/bin/env bash
# =========================================================
# embryo_gained_lost_heatmaps.sh
# ---------------------------------------------------------
# Lifts GSE101571 8C and ICM ATAC-seq data (hg19) over to hg38,
# then generates heatmaps of embryo ATAC signal over the
# gained/lost peaks of every condition vs DMSO.
# =========================================================

set -euo pipefail

BASE_DIR="path/to/ATAC_analysis"
PROC_DIR="${BASE_DIR}/processed_multiNorm"
DIFF_DIR="${PROC_DIR}/diff_atac"
MATRIX_DIR="${PROC_DIR}/matrices"
HEATMAP_DIR="${PROC_DIR}/heatmaps_diff"
GENOME_DIR="path/to/hg38"
CHROM_SIZES="${GENOME_DIR}/hg38.chrom.sizes"
CHAIN="${GENOME_DIR}/hg19ToHg38.over.chain.gz"
SAMPLES="${BASE_DIR}/samples.tsv"
THREADS=6

EMBRYO_DIR="${BASE_DIR}/embryo_data"
mkdir -p "${EMBRYO_DIR}" "${HEATMAP_DIR}"

WIG_8C="${EMBRYO_DIR}/GSE101571_8cell_2pn_100bp_rpkm.wig.gz"
WIG_ICM="${EMBRYO_DIR}/GSE101571_icm_2pn_100bp_rpkm.wig.gz"
BW_8C="${EMBRYO_DIR}/GSE101571_8cell_rpkm_hg38.bw"
BW_ICM="${EMBRYO_DIR}/GSE101571_icm_rpkm_hg38.bw"

# --- Chain file -------------------------------------------------------------
if [[ ! -f "${CHAIN}" ]]; then
    echo "Downloading chain file..."
    wget -q https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz \
      -O "${CHAIN}"
fi

# --- Step A: wig.gz -> liftover hg19->hg38 -> bigWig -------------------------
echo "=== A. wig.gz -> liftover hg19->hg38 -> bigWig ==="

convert_and_liftover () {
    local WIG_GZ=$1
    local BW_OUT=$2
    local LABEL=$3
    local TMP="${EMBRYO_DIR}/${LABEL}_tmp"

    if [[ -f "${BW_OUT}" ]]; then
        echo "  ${LABEL}: hg38 bigWig exists - skipping."
        return
    fi

    local HG19_SIZES="${EMBRYO_DIR}/hg19.chrom.sizes"
    if [[ ! -f "${HG19_SIZES}" ]]; then
        echo "  Downloading hg19 chrom sizes..."
        wget -q https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes \
          -O "${HG19_SIZES}"
    fi

    echo "  ${LABEL}: wig.gz -> bigWig (hg19)..."
    gunzip -c "${WIG_GZ}" > "${TMP}.hg19.wig"
    wigToBigWig "${TMP}.hg19.wig" "${HG19_SIZES}" "${TMP}.hg19.bw"
    rm -f "${TMP}.hg19.wig"

    echo "  ${LABEL}: bigWig -> bedGraph..."
    bigWigToBedGraph "${TMP}.hg19.bw" "${TMP}.hg19.bedGraph"
    rm -f "${TMP}.hg19.bw"

    echo "  ${LABEL}: liftover hg19 -> hg38..."
    liftOver \
      "${TMP}.hg19.bedGraph" \
      "${CHAIN}" \
      "${TMP}.hg38.bedGraph" \
      "${TMP}.unmapped.bed"

    UNMAPPED=$(wc -l < "${TMP}.unmapped.bed" || echo 0)
    echo "  ${LABEL}: unmapped lines: ${UNMAPPED}"

    echo "  ${LABEL}: sort + clean..."
    sort -k1,1 -k2,2n "${TMP}.hg38.bedGraph" \
      | awk 'BEGIN{OFS="\t"} $1 ~ /^chr([0-9]+|X|Y)$/ {print}' \
      | bedtools merge -i - -c 4 -o mean \
      > "${TMP}.hg38.clean.bedGraph"

    echo "  ${LABEL}: bedGraph -> bigWig (hg38)..."
    bedGraphToBigWig "${TMP}.hg38.clean.bedGraph" "${CHROM_SIZES}" "${BW_OUT}"

    rm -f "${TMP}.hg19.bedGraph" "${TMP}.hg38.bedGraph" \
          "${TMP}.hg38.clean.bedGraph" "${TMP}.unmapped.bed"

    echo "  ${LABEL}: done -> ${BW_OUT}"
}

convert_and_liftover "${WIG_8C}"  "${BW_8C}"  "8cell"
convert_and_liftover "${WIG_ICM}" "${BW_ICM}" "ICM"

# --- Step B: heatmap for all conditions vs DMSO -----------------------------
echo "=== B. Heatmap generation - all conditions vs DMSO ==="

CONDITIONS=$(tail -n +2 "${SAMPLES}" | awk -F'\t' '$2 != "DMSO" {print $2}' | sort -u)

for COND in ${CONDITIONS}; do
    for DIRECTION in gained lost; do

        BED="${DIFF_DIR}/${COND}_vs_DMSO_${DIRECTION}.bed"

        if [[ ! -s "${BED}" ]]; then
            echo "  ${COND}_vs_DMSO_${DIRECTION}: empty BED - skipping."
            continue
        fi

        N_PEAKS=$(wc -l < "${BED}")
        MAT="${MATRIX_DIR}/${COND}_vs_DMSO_${DIRECTION}_embryo.gz"
        PDF="${HEATMAP_DIR}/${COND}_vs_DMSO_${DIRECTION}_embryo_heatmap.pdf"

        if [[ -f "${PDF}" ]]; then
            echo "  ${COND}_vs_DMSO_${DIRECTION}: heatmap exists - skipping."
            continue
        fi

        if [[ ! -f "${MAT}" ]]; then
            echo "  Computing matrix: ${COND}_vs_DMSO_${DIRECTION} (n=${N_PEAKS})"
            computeMatrix reference-point \
              --referencePoint center \
              -b 2000 -a 2000 \
              -R "${BED}" \
              -S "${BW_8C}" "${BW_ICM}" \
              --skipZeros \
              --missingDataAsZero \
              -p ${THREADS} \
              -o "${MAT}"
        fi

        if [[ "${DIRECTION}" == "gained" ]]; then
            COLORMAP="Reds"
            TITLE="${COND} vs DMSO - gained (n=${N_PEAKS})"
        else
            COLORMAP="Blues"
            TITLE="${COND} vs DMSO - lost (n=${N_PEAKS})"
        fi

        echo "  Plotting: ${COND}_vs_DMSO_${DIRECTION}"
        plotHeatmap -m "${MAT}" \
          --colorMap ${COLORMAP} \
          --whatToShow 'heatmap and colorbar' \
          --samplesLabel "8C" "ICM" \
          --refPointLabel "center" \
          --xAxisLabel "Distance from center (bp)" \
          --zMin 0 \
          --sortUsing mean \
          --sortRegions descend \
          --missingDataColor 'white' \
          --heatmapHeight 18 \
          --heatmapWidth 3 \
          --plotTitle "${TITLE}" \
          -out "${PDF}"

        echo "  Saved: ${PDF}"

    done
done

echo ""
echo "Done."
echo "  hg38 bigWigs : ${EMBRYO_DIR}/*_hg38.bw"
echo "  Heatmaps     : ${HEATMAP_DIR}/*_embryo_heatmap.pdf"
