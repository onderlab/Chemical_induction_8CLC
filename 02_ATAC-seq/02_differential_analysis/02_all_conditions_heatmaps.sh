#!/usr/bin/env bash
# =========================================================
# all_conditions_heatmaps.sh
# ---------------------------------------------------------
# Each condition vs DMSO: gained and lost ATAC peak heatmaps.
# Uses hkTSS-normalized bigWigs. For each comparison, only DMSO
# plus the relevant condition's bigWigs are shown.
# =========================================================

set -euo pipefail

BASE_DIR="path/to/ATAC_analysis"
PROC_DIR="${BASE_DIR}/processed_multiNorm"
COUNT_DIR="${PROC_DIR}/counts"
BW_DIR="${PROC_DIR}/bigwig"
MATRIX_DIR="${PROC_DIR}/matrices"
DIFF_DIR="${PROC_DIR}/diff_atac"
HEATMAP_DIR="${PROC_DIR}/heatmaps_diff"
SCRIPT_DIR="${BASE_DIR}/scripts"
SAMPLES="${BASE_DIR}/samples.tsv"
THREADS=6

mkdir -p "${DIFF_DIR}" "${HEATMAP_DIR}"

# --- Step A: DESeq2 (if not run yet) ----------------------------------------
echo "=== A. Differential ATAC: all conditions vs DMSO ==="

if [[ -f "${DIFF_DIR}/diff_atac_summary.csv" ]]; then
    echo "  DESeq2 results exist - skipping."
else
    Rscript "${SCRIPT_DIR}/diff_atac_vs_dmso.R" \
      "${COUNT_DIR}/consensus_peak_counts.txt" \
      "${SAMPLES}" \
      "${COUNT_DIR}/hkTSS_scale_factors.tsv" \
      "${DIFF_DIR}"
fi

# --- Step B: find DMSO sample names -----------------------------------------
DMSO_SAMPLES=()
while IFS=$'\t' read -r SAMPLE CONDITION REP BAM; do
    [[ "${CONDITION}" == "DMSO" ]] && DMSO_SAMPLES+=("${SAMPLE}")
done < <(tail -n +2 "${SAMPLES}")

echo "  DMSO samples: ${DMSO_SAMPLES[*]}"

# --- Step C: loop over each condition ---------------------------------------
echo "=== B. Heatmap generation ==="

# Unique non-DMSO conditions
CONDITIONS=$(tail -n +2 "${SAMPLES}" | awk -F'\t' '$2 != "DMSO" {print $2}' | sort -u)

for COND in ${CONDITIONS}; do

    echo ""
    echo "-- ${COND} vs DMSO --"

    # Sample names for this condition
    COND_SAMPLES=()
    while IFS=$'\t' read -r SAMPLE CONDITION REP BAM; do
        [[ "${CONDITION}" == "${COND}" ]] && COND_SAMPLES+=("${SAMPLE}")
    done < <(tail -n +2 "${SAMPLES}")

    # bigWig list: DMSO first, then the condition
    BWFILES=()
    LABELS=()
    for S in "${DMSO_SAMPLES[@]}" "${COND_SAMPLES[@]}"; do
        BW="${BW_DIR}/${S}.hkTSS.bw"
        if [[ ! -f "${BW}" ]]; then
            echo "  ERROR: ${BW} not found - complete the pipeline first."
            exit 1
        fi
        BWFILES+=("${BW}")
        LABELS+=("${S}")
    done
    LABEL_STR="${LABELS[*]}"

    for DIRECTION in gained lost; do

        BED="${DIFF_DIR}/${COND}_vs_DMSO_${DIRECTION}.bed"
        MAT="${MATRIX_DIR}/${COND}_vs_DMSO_${DIRECTION}_hkTSS.gz"
        PDF="${HEATMAP_DIR}/${COND}_vs_DMSO_${DIRECTION}_hkTSS_heatmap.pdf"

        if [[ ! -s "${BED}" ]]; then
            echo "  ${COND}_vs_DMSO_${DIRECTION}: empty BED - skipping."
            continue
        fi

        N_PEAKS=$(wc -l < "${BED}")
        echo "  ${DIRECTION}: ${N_PEAKS} peaks"

        if [[ ! -f "${MAT}" ]]; then
            echo "  Computing matrix: ${COND}_vs_DMSO_${DIRECTION}"
            computeMatrix reference-point \
              --referencePoint center \
              -b 2000 -a 2000 \
              -R "${BED}" \
              -S "${BWFILES[@]}" \
              --skipZeros \
              --missingDataAsZero \
              -p ${THREADS} \
              -o "${MAT}"
        fi

        if [[ "${DIRECTION}" == "gained" ]]; then
            COLORMAP="Reds"
            TITLE="ATAC signal gained (n=${N_PEAKS})"
        else
            COLORMAP="Blues"
            TITLE="ATAC signal lost (n=${N_PEAKS})"
        fi

        if [[ -f "${PDF}" ]]; then
            echo "  ${COND}_vs_DMSO_${DIRECTION}: heatmap exists - skipping."
            continue
        fi

        echo "  Plotting: ${COND}_vs_DMSO_${DIRECTION}"
        plotHeatmap -m "${MAT}" \
          --colorMap ${COLORMAP} \
          --whatToShow 'heatmap and colorbar' \
          --samplesLabel ${LABEL_STR} \
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
echo "  Heatmaps: ${HEATMAP_DIR}/"
