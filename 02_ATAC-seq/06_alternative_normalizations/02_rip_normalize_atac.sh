#!/usr/bin/env bash
# =========================================================
# rip_normalize_atac.sh
# ---------------------------------------------------------
# ATAC-seq Reads-in-Peaks (RIP) normalization.
# Outputs: rip_normalization_factors.tsv, *.RIP.bw, ZGA heatmap.
# =========================================================

set -euo pipefail

BASE_DIR="path/to/ATAC_analysis"
PROC_DIR="${BASE_DIR}/processed_multiNorm"
SHIFT_DIR="${PROC_DIR}/shifted_bam"
COUNT_DIR="${PROC_DIR}/counts"
BW_DIR="${PROC_DIR}/bigwig"
MATRIX_DIR="${PROC_DIR}/matrices"
HEATMAP_DIR="${PROC_DIR}/heatmaps_diff"
BED_DIR="${PROC_DIR}/beds"
SAMPLES="${BASE_DIR}/samples.tsv"
THREADS=6
BINSIZE=10
EFFECTIVE_GENOME_SIZE=2913022398

RIP_DIR="${PROC_DIR}/rip_norm"
mkdir -p "${RIP_DIR}/tmp" "${BW_DIR}" "${MATRIX_DIR}" "${HEATMAP_DIR}"

PEAKS="${COUNT_DIR}/consensus_peaks.bed"
COUNTS_TSV="${RIP_DIR}/rip_counts.tsv"
NORM_TSV="${RIP_DIR}/rip_normalization_factors.tsv"

# --- Step 1: count fragments in peaks ---------------------------------------
echo "=== 1. Reads-in-Peaks counting ==="

if [[ -f "${NORM_TSV}" ]]; then
    echo "  RIP normalization factors exist - skipping."
else
    SORTED_PEAKS="${RIP_DIR}/tmp/consensus_peaks.sorted.bed"
    sort -k1,1 -k2,2n "${PEAKS}" > "${SORTED_PEAKS}"

    echo -e "sample\tbam\tfragments_in_peaks" > "${COUNTS_TSV}"

    tail -n +2 "${SAMPLES}" | while IFS=$'\t' read -r SAMPLE CONDITION REP BAM; do
        BAM_SHIFT="${SHIFT_DIR}/${SAMPLE}.shifted.bam"

        if [[ ! -f "${BAM_SHIFT}" ]]; then
            echo "  ERROR: ${BAM_SHIFT} not found"
            exit 1
        fi

        echo "  Counting: ${SAMPLE}"

        # Name-sort BAM for bedpe conversion
        NAMESORTED="${RIP_DIR}/tmp/${SAMPLE}.namesorted.bam"
        if [[ ! -f "${NAMESORTED}" ]]; then
            samtools sort -n -@ ${THREADS} \
              -o "${NAMESORTED}" \
              "${SHIFT_DIR}/${SAMPLE}.shifted.bam"
        fi

        bedtools bamtobed -bedpe -i "${NAMESORTED}" \
            | awk '$1==$4 {print $1, $2, $6}' OFS='\t' \
            | bedtools intersect -u -a - -b "${SORTED_PEAKS}" \
            | wc -l > "${RIP_DIR}/tmp/${SAMPLE}.count"

        COUNT=$(cat "${RIP_DIR}/tmp/${SAMPLE}.count" | tr -d ' ')
        echo -e "${SAMPLE}\t${BAM_SHIFT}\t${COUNT}" >> "${COUNTS_TSV}"
        echo "  ${SAMPLE}: ${COUNT} fragments in peaks"
        rm -f "${NAMESORTED}"
    done

    # Compute scale factors
    python3 - <<EOF
import pandas as pd

df = pd.read_csv("${COUNTS_TSV}", sep="\t")
median = df["fragments_in_peaks"].median()

df["scale_per_million"] = 1_000_000 / df["fragments_in_peaks"]
df["scale_to_median"]   = median / df["fragments_in_peaks"]

print("\nRIP normalization factors:")
print(df[["sample","fragments_in_peaks","scale_to_median"]].to_string(index=False))

df.to_csv("${NORM_TSV}", sep="\t", index=False)
EOF

    echo "  Norm factors: ${NORM_TSV}"
fi

# --- Step 2: produce RIP-normalized bigWigs ---------------------------------
echo ""
echo "=== 2. RIP-normalized bigWig generation ==="

# Read the scale_to_median column
SKIP_BW=true
tail -n +2 "${NORM_TSV}" | while IFS=$'\t' read -r SAMPLE BAM FIP SPM STM; do
    [[ ! -f "${BW_DIR}/${SAMPLE}.RIP.bw" ]] && exit 1
    true
done && SKIP_BW=true || SKIP_BW=false

if ${SKIP_BW}; then
    echo "  All RIP bigWigs exist - skipping."
else
    tail -n +2 "${NORM_TSV}" | while IFS=$'\t' read -r SAMPLE BAM FIP SPM STM; do
        OUT="${BW_DIR}/${SAMPLE}.RIP.bw"

        if [[ -f "${OUT}" ]]; then
            echo "  ${SAMPLE} RIP bigWig exists - skipping."
            continue
        fi

        echo "  bigWig: ${SAMPLE} (scale_to_median=${STM})"
        bamCoverage \
            -b "${SHIFT_DIR}/${SAMPLE}.shifted.bam" \
            -o "${OUT}" \
            --binSize ${BINSIZE} \
            --normalizeUsing None \
            --scaleFactor "${STM}" \
            --extendReads \
            -p ${THREADS}
    done
fi

# --- Step 3: ZGA TSS heatmap (RIP-normalized) -------------------------------
echo ""
echo "=== 3. ZGA TSS heatmap (RIP normalized) ==="

ZGA_BED="${BED_DIR}/zga_tss_1bp.bed"
MAT="${MATRIX_DIR}/zgaTSS_RIP.gz"
PDF="${HEATMAP_DIR}/zgaTSS_RIP_heatmap.pdf"
PROFILE_PDF="${HEATMAP_DIR}/zgaTSS_RIP_profile.pdf"
PROFILE_TSV="${HEATMAP_DIR}/zgaTSS_RIP_profile.tsv"

# bigWig list (read directly to avoid subshell scope issues)
BWFILES=()
LABELS=()
while IFS=$'\t' read -r SAMPLE BAM FIP SPM STM; do
    BWFILES+=("${BW_DIR}/${SAMPLE}.RIP.bw")
    LABELS+=("${SAMPLE}")
done < <(tail -n +2 "${NORM_TSV}")
LABEL_STR="${LABELS[*]}"

if [[ ! -f "${MAT}" ]]; then
    N_ZGA=$(wc -l < "${ZGA_BED}")
    echo "  Computing matrix: ZGA TSS (n=${N_ZGA})"
    computeMatrix reference-point \
      --referencePoint TSS \
      -b 2000 -a 2000 \
      -R "${ZGA_BED}" \
      -S "${BWFILES[@]}" \
      --skipZeros \
      --missingDataAsZero \
      -p ${THREADS} \
      -o "${MAT}"
fi

if [[ ! -f "${PROFILE_PDF}" ]]; then
    echo "  plotProfile: ZGA TSS RIP"
    plotProfile \
      -m "${MAT}" \
      --perGroup \
      --samplesLabel ${LABEL_STR} \
      --refPointLabel "TSS" \
      --outFileNameData "${PROFILE_TSV}" \
      -out "${PROFILE_PDF}"
fi

if [[ ! -f "${PDF}" ]]; then
    N_ZGA=$(wc -l < "${ZGA_BED}")
    echo "  plotHeatmap: ZGA TSS RIP"
    plotHeatmap -m "${MAT}" \
      --colorMap Reds \
      --whatToShow 'heatmap and colorbar' \
      --samplesLabel ${LABEL_STR} \
      --refPointLabel "TSS" \
      --xAxisLabel "Distance from TSS (bp)" \
      --zMin 0 \
      --sortUsing mean \
      --sortRegions descend \
      --missingDataColor 'white' \
      --heatmapHeight 18 \
      --heatmapWidth 3 \
      --plotTitle "ZGA TSS - RIP normalization (n=${N_ZGA})" \
      -out "${PDF}"
    echo "  Saved: ${PDF}"
fi

echo ""
echo "Done."
echo "  Norm factors : ${NORM_TSV}"
echo "  bigWigs      : ${BW_DIR}/*.RIP.bw"
echo "  ZGA profile  : ${PROFILE_PDF}"
echo "  ZGA heatmap  : ${PDF}"
