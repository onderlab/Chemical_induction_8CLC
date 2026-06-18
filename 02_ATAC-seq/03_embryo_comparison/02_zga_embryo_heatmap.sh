#!/usr/bin/env bash
# =========================================================
# zga_embryo_heatmap.sh
# ---------------------------------------------------------
# Shows 8C and ICM embryo ATAC-seq signal over ZGA gene TSS.
# Requires the hg38 embryo bigWigs produced by
# 01_embryo_liftover_gained_lost_heatmaps.sh.
# =========================================================

set -euo pipefail

BASE_DIR="path/to/ATAC_analysis"
PROC_DIR="${BASE_DIR}/processed_multiNorm"
BED_DIR="${PROC_DIR}/beds"
MATRIX_DIR="${PROC_DIR}/matrices"
HEATMAP_DIR="${PROC_DIR}/heatmaps_diff"
EMBRYO_DIR="${BASE_DIR}/embryo_data"
GENOME_DIR="path/to/hg38"
CHROM_SIZES="${GENOME_DIR}/hg38.chrom.sizes"
CHAIN="${GENOME_DIR}/hg19ToHg38.over.chain.gz"
THREADS=6

mkdir -p "${HEATMAP_DIR}"

BW_8C="${EMBRYO_DIR}/GSE101571_8cell_rpkm_hg38.bw"
BW_ICM="${EMBRYO_DIR}/GSE101571_icm_rpkm_hg38.bw"
ZGA_BED="${BED_DIR}/zga_tss_1bp.bed"

# --- Check inputs -----------------------------------------------------------
for F in "${BW_8C}" "${BW_ICM}" "${ZGA_BED}"; do
    if [[ ! -f "${F}" ]]; then
        echo "ERROR: file not found: ${F}"
        echo "Run 01_embryo_liftover_gained_lost_heatmaps.sh first."
        exit 1
    fi
done

N_ZGA=$(wc -l < "${ZGA_BED}")
echo "ZGA gene count: ${N_ZGA}"

# --- computeMatrix ----------------------------------------------------------
MAT="${MATRIX_DIR}/zgaTSS_embryo.gz"

if [[ ! -f "${MAT}" ]]; then
    echo "=== computeMatrix: ZGA TSS +/- 2kb ==="
    computeMatrix reference-point \
      --referencePoint TSS \
      -b 2000 -a 2000 \
      -R "${ZGA_BED}" \
      -S "${BW_8C}" "${BW_ICM}" \
      --skipZeros \
      --missingDataAsZero \
      -p ${THREADS} \
      -o "${MAT}"
else
    echo "Matrix exists - skipping computeMatrix."
fi

# --- plotProfile ------------------------------------------------------------
PROFILE_PDF="${HEATMAP_DIR}/zgaTSS_embryo_profile.pdf"
PROFILE_TSV="${HEATMAP_DIR}/zgaTSS_embryo_profile.tsv"

if [[ ! -f "${PROFILE_PDF}" ]]; then
    echo "=== plotProfile: ZGA TSS ==="
    plotProfile \
      -m "${MAT}" \
      --perGroup \
      --samplesLabel "8C" "ICM" \
      --refPointLabel "TSS" \
      --outFileNameData "${PROFILE_TSV}" \
      -out "${PROFILE_PDF}"
fi

# --- plotHeatmap ------------------------------------------------------------
HEATMAP_PDF="${HEATMAP_DIR}/zgaTSS_embryo_heatmap.pdf"

if [[ ! -f "${HEATMAP_PDF}" ]]; then
    echo "=== plotHeatmap: ZGA TSS ==="
    plotHeatmap -m "${MAT}" \
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
      --plotTitle "ZGA gene TSS (n=${N_ZGA})" \
      -out "${HEATMAP_PDF}"
fi

echo ""
echo "Done."
echo "  Profile : ${PROFILE_PDF}"
echo "  Heatmap : ${HEATMAP_PDF}"
