#!/usr/bin/env bash
# =========================================================
# run_downsampled_pipeline.sh
# ---------------------------------------------------------
# Downsamples all shifted BAMs to the sample with the fewest
# reads, then re-runs CPM/RPGC/stablePeak/hkTSS normalizations
# and heatmaps. Used as a control to confirm that observed
# accessibility differences are not driven by library depth.
#
# NOTE: TARGET_READS and the TOTAL_READS table below are
# specific to this dataset. Replace with your own flagstat
# counts before running.
# =========================================================

set -euo pipefail

BASE_DIR="path/to/ATAC_analysis"
PROC_DIR="${BASE_DIR}/processed_multiNorm"
SHIFT_DIR="${PROC_DIR}/shifted_bam"
GENOME_DIR="path/to/hg38"
CHROM_SIZES="${GENOME_DIR}/hg38.chrom.sizes"
SAMPLES="${BASE_DIR}/samples.tsv"
SCRIPT_DIR="${BASE_DIR}/scripts"
THREADS=6
EFFECTIVE_GENOME_SIZE=2913022398

# Downsampled output directories - kept separate from the originals
DS_DIR="${PROC_DIR}/downsampled"
DS_SHIFT_DIR="${DS_DIR}/shifted_bam"
DS_COUNT_DIR="${DS_DIR}/counts"
DS_BW_DIR="${DS_DIR}/bigwig"
DS_MATRIX_DIR="${DS_DIR}/matrices"
DS_HEATMAP_DIR="${DS_DIR}/heatmaps"
DS_QC_DIR="${DS_DIR}/qc"
DS_BED_DIR="${PROC_DIR}/beds"   # BED files unchanged - reuse

mkdir -p "${DS_SHIFT_DIR}" "${DS_COUNT_DIR}" "${DS_BW_DIR}" \
         "${DS_MATRIX_DIR}" "${DS_HEATMAP_DIR}" "${DS_QC_DIR}"

# --- Downsampling parameters (dataset-specific) -----------------------------
TARGET_READS=31857818   # lowest-read sample

declare -A TOTAL_READS
TOTAL_READS["DMSO_Rep1"]=80602728
TOTAL_READS["DMSO_Rep2"]=76076240
TOTAL_READS["DUX4sorted_Rep1"]=136670879
TOTAL_READS["DUX4sorted_Rep2"]=92956670
TOTAL_READS["PRBJN_Rep1"]=64061521
TOTAL_READS["PRBJN_Rep2"]=47480446
TOTAL_READS["PRBJNsorted_Rep1"]=43282920
TOTAL_READS["PRBJNsorted_Rep2"]=31857818
TOTAL_READS["PY60_Rep1"]=41545720
TOTAL_READS["PY60_Rep2"]=46781384

# --- Step 1: Downsampling ---------------------------------------------------
echo "=== 1. Downsampling -> ${TARGET_READS} reads ==="

tail -n +2 "${SAMPLES}" | while IFS=$'\t' read -r SAMPLE CONDITION REP BAM; do
    IN_BAM="${SHIFT_DIR}/${SAMPLE}.shifted.bam"
    OUT_BAM="${DS_SHIFT_DIR}/${SAMPLE}.shifted.bam"

    if [[ -f "${OUT_BAM}" && -f "${OUT_BAM}.bai" ]]; then
        echo "  ${SAMPLE}: exists - skipping."
        continue
    fi

    TOTAL=${TOTAL_READS[$SAMPLE]}

    if [[ "${TOTAL}" -le "${TARGET_READS}" ]]; then
        # Already below the target - copy as is
        echo "  ${SAMPLE}: ${TOTAL} <= ${TARGET_READS} - copying."
        cp "${IN_BAM}" "${OUT_BAM}"
        cp "${IN_BAM}.bai" "${OUT_BAM}.bai"
    else
        # Compute downsampling fraction
        # samtools -s: integer part = seed, decimal part = fraction to keep
        FRAC_STR=$(python3 -c "t=${TARGET_READS};tot=${TOTAL_READS[$SAMPLE]};print(f'{t/tot:.4f}')")
        FRAC_DEC="${FRAC_STR#*.}"
        echo "  ${SAMPLE}: ${TOTAL} -> ${TARGET_READS} (frac=${FRAC_STR})"
        samtools view -@ ${THREADS} -b \
          -s "42.${FRAC_DEC}" \
          -o "${OUT_BAM}" \
          "${IN_BAM}"
        samtools index -@ ${THREADS} "${OUT_BAM}"
    fi

    samtools flagstat "${OUT_BAM}" > "${DS_QC_DIR}/${SAMPLE}.flagstat.txt"
done

echo ""
echo "=== Downsampled read counts ==="
for f in "${DS_QC_DIR}"/*.flagstat.txt; do
    echo -n "  $(basename $f .flagstat.txt): "
    grep "mapped (" "$f" | awk '{print $1}'
done

# --- Step 2: Peak calling (consensus peaks unchanged - skip) ----------------
CONSENSUS_BED="${PROC_DIR}/counts/consensus_peaks.bed"
CONSENSUS_SAF="${PROC_DIR}/counts/consensus_peaks.saf"

if [[ ! -f "${CONSENSUS_BED}" ]]; then
    echo "ERROR: consensus peaks not found: ${CONSENSUS_BED}"
    exit 1
fi
echo ""
echo "=== 2. Consensus peaks exist - skipping ==="

# --- Step 3: featureCounts --------------------------------------------------
echo ""
echo "=== 3. featureCounts (downsampled) ==="

COUNTS_FILE="${DS_COUNT_DIR}/consensus_peak_counts.txt"

if [[ -f "${COUNTS_FILE}" ]]; then
    echo "  Count file exists - skipping."
else
    BAMS=()
    while IFS=$'\t' read -r SAMPLE CONDITION REP BAM; do
        [[ "${SAMPLE}" == "sample" ]] && continue
        BAMS+=("${DS_SHIFT_DIR}/${SAMPLE}.shifted.bam")
    done < "${SAMPLES}"

    featureCounts \
      -a "${CONSENSUS_SAF}" \
      -F SAF \
      -T ${THREADS} \
      -p -B -C \
      -o "${COUNTS_FILE}" \
      "${BAMS[@]}"
fi

# --- Step 4: R normalization (stablePeak + HK) ------------------------------
echo ""
echo "=== 4. R normalization ==="

if [[ -f "${DS_COUNT_DIR}/stablePeak_scale_factors.tsv" && \
      -f "${DS_QC_DIR}/PCA_stablePeakNorm.pdf" ]]; then
    echo "  Scale factors exist - skipping."
else
    Rscript "${SCRIPT_DIR}/atac_multinorm.R" \
      "${COUNTS_FILE}" \
      "${SAMPLES}" \
      "${DS_BED_DIR}/hk_tss_1bp.bed" \
      "${DS_COUNT_DIR}" \
      "${DS_QC_DIR}"
fi

# --- Step 5: bigWig - CPM + RPGC --------------------------------------------
echo ""
echo "=== 5. CPM + RPGC bigWig ==="

tail -n +2 "${SAMPLES}" | while IFS=$'\t' read -r SAMPLE CONDITION REP BAM; do
    DS_BAM="${DS_SHIFT_DIR}/${SAMPLE}.shifted.bam"

    if [[ ! -f "${DS_BW_DIR}/${SAMPLE}.CPM.bw" ]]; then
        echo "  CPM: ${SAMPLE}"
        bamCoverage -b "${DS_BAM}" -o "${DS_BW_DIR}/${SAMPLE}.CPM.bw" \
          --binSize 10 --normalizeUsing CPM --extendReads -p ${THREADS}
    fi

    if [[ ! -f "${DS_BW_DIR}/${SAMPLE}.RPGC.bw" ]]; then
        echo "  RPGC: ${SAMPLE}"
        bamCoverage -b "${DS_BAM}" -o "${DS_BW_DIR}/${SAMPLE}.RPGC.bw" \
          --binSize 10 --normalizeUsing RPGC \
          --effectiveGenomeSize ${EFFECTIVE_GENOME_SIZE} --extendReads -p ${THREADS}
    fi
done

# --- Step 6: bigWig - stablePeak --------------------------------------------
echo ""
echo "=== 6. stablePeak bigWig ==="

tail -n +2 "${DS_COUNT_DIR}/stablePeak_scale_factors.tsv" | \
while IFS=$'\t' read -r SAMPLE SIZEFACTOR SCALEFACTOR; do
    DS_BAM="${DS_SHIFT_DIR}/${SAMPLE}.shifted.bam"
    OUT="${DS_BW_DIR}/${SAMPLE}.stablePeak.bw"
    if [[ ! -f "${OUT}" ]]; then
        echo "  stablePeak: ${SAMPLE} (scale=${SCALEFACTOR})"
        bamCoverage -b "${DS_BAM}" -o "${OUT}" \
          --binSize 10 --normalizeUsing None \
          --scaleFactor "${SCALEFACTOR}" --extendReads -p ${THREADS}
    fi
done

# --- Step 7: hkTSS scale factors + bigWig -----------------------------------
echo ""
echo "=== 7. hkTSS normalization ==="

if [[ ! -f "${DS_COUNT_DIR}/hkTSS_scale_factors.tsv" ]]; then
    Rscript "${SCRIPT_DIR}/hkTSS_normalize.R" \
      "${DS_BED_DIR}/hk_tss_1bp.bed" \
      "${DS_SHIFT_DIR}" \
      "${SAMPLES}" \
      "${DS_COUNT_DIR}" \
      "${DS_QC_DIR}"
fi

tail -n +2 "${DS_COUNT_DIR}/hkTSS_scale_factors.tsv" | \
while IFS=$'\t' read -r SAMPLE _READS SCALEFACTOR; do
    DS_BAM="${DS_SHIFT_DIR}/${SAMPLE}.shifted.bam"
    OUT="${DS_BW_DIR}/${SAMPLE}.hkTSS.bw"
    if [[ ! -f "${OUT}" ]]; then
        echo "  hkTSS: ${SAMPLE} (scale=${SCALEFACTOR})"
        bamCoverage -b "${DS_BAM}" -o "${OUT}" \
          --binSize 10 --normalizeUsing None \
          --scaleFactor "${SCALEFACTOR}" --extendReads -p ${THREADS}
    fi
done

# --- Step 8: computeMatrix + plotHeatmap (ZGA TSS) --------------------------
echo ""
echo "=== 8. ZGA TSS heatmaps ==="

ZGA_BED="${DS_BED_DIR}/zga_tss_1bp.bed"
HK_BED="${DS_BED_DIR}/hk_tss_1bp.bed"
ALL_BED="${DS_BED_DIR}/all_gene_tss_1bp.bed"

SAMPLE_LABELS=$(cut -f1 "${SAMPLES}" | tail -n +2 | tr '\n' ' ')

for MODE in CPM RPGC stablePeak hkTSS; do

    BWFILES=()
    if [[ "${MODE}" == "stablePeak" ]]; then
        while IFS=$'\t' read -r S _SF _SC; do
            BWFILES+=("${DS_BW_DIR}/${S}.${MODE}.bw")
        done < <(tail -n +2 "${DS_COUNT_DIR}/stablePeak_scale_factors.tsv")
    elif [[ "${MODE}" == "hkTSS" ]]; then
        while IFS=$'\t' read -r S _R _SC; do
            BWFILES+=("${DS_BW_DIR}/${S}.hkTSS.bw")
        done < <(tail -n +2 "${DS_COUNT_DIR}/hkTSS_scale_factors.tsv")
    else
        while IFS=$'\t' read -r SAMPLE CONDITION REP BAM; do
            [[ "${SAMPLE}" == "sample" ]] && continue
            BWFILES+=("${DS_BW_DIR}/${SAMPLE}.${MODE}.bw")
        done < "${SAMPLES}"
    fi

    for REG in allTSS hkTSS zgaTSS; do
        case "${REG}" in
            allTSS) BED="${ALL_BED}" ;;
            hkTSS)  BED="${HK_BED}" ;;
            zgaTSS) BED="${ZGA_BED}" ;;
        esac

        MAT="${DS_MATRIX_DIR}/${REG}_${MODE}.gz"
        PDF="${DS_HEATMAP_DIR}/${REG}_${MODE}_heatmap.pdf"
        PROFILE_PDF="${DS_HEATMAP_DIR}/${REG}_${MODE}_profile.pdf"
        PROFILE_TSV="${DS_MATRIX_DIR}/${REG}_${MODE}_profile.tsv"

        if [[ ! -f "${MAT}" ]]; then
            echo "  Matrix: ${REG}_${MODE}"
            computeMatrix reference-point \
              --referencePoint TSS -b 2000 -a 2000 \
              -R "${BED}" -S "${BWFILES[@]}" \
              --skipZeros --missingDataAsZero -p ${THREADS} \
              -o "${MAT}"
        fi

        if [[ ! -f "${PROFILE_PDF}" ]]; then
            plotProfile -m "${MAT}" --perGroup \
              --samplesLabel ${SAMPLE_LABELS} \
              --refPointLabel "TSS" \
              --outFileNameData "${PROFILE_TSV}" \
              -out "${PROFILE_PDF}"
        fi

        if [[ "${REG}" == "zgaTSS" && ! -f "${PDF}" ]]; then
            plotHeatmap -m "${MAT}" \
              --colorMap Reds \
              --whatToShow 'heatmap and colorbar' \
              --samplesLabel ${SAMPLE_LABELS} \
              --refPointLabel "TSS" \
              --xAxisLabel "Distance from TSS (bp)" \
              --zMin 0 --sortUsing mean --sortRegions descend \
              --missingDataColor 'white' \
              --heatmapHeight 18 --heatmapWidth 3 \
              --plotTitle "ZGA TSS - ${MODE} (downsampled)" \
              -out "${PDF}"
            echo "  Saved: ${PDF}"
        fi
    done
done

echo ""
echo "Done."
echo "  Downsampled BAMs : ${DS_SHIFT_DIR}/"
echo "  bigWigs          : ${DS_BW_DIR}/"
echo "  Heatmaps         : ${DS_HEATMAP_DIR}/"
echo "  QC               : ${DS_QC_DIR}/"
