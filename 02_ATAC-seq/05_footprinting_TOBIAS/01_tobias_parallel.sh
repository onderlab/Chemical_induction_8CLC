#!/usr/bin/env bash
# =========================================================
# tobias_parallel.sh
# ---------------------------------------------------------
# TOBIAS digital footprinting, run in parallel per sample.
#   1. ATACorrect  - correct Tn5 insertion bias
#   2. ScoreBigwig - compute footprint scores
#   3. BINDetect   - differential TF binding (condition vs DMSO)
# =========================================================

export OBJC_DISABLE_INITIALIZE_FORK_SAFETY=YES
export PYTHONPATH=""

BASE_DIR="path/to/ATAC_analysis"
PROC_DIR="${BASE_DIR}/processed_multiNorm"
SHIFT_DIR="${PROC_DIR}/shifted_bam"
TOBIAS_DIR="${PROC_DIR}/tobias"
GENOME_FA="path/to/hg38/hg38.fa"
BLACKLIST="path/to/hg38/hg38-blacklist.v2.bed"
CONSENSUS_BED="${PROC_DIR}/counts/consensus_peaks.bed"
SAMPLES="${BASE_DIR}/samples.tsv"
MOTIF_DIR="${PROC_DIR}/motif_analysis"
JASPAR_MEME="${MOTIF_DIR}/JASPAR2024_vertebrates.meme"

mkdir -p "${TOBIAS_DIR}"

# Sample list
SAMPLE_LIST=()
while IFS=$'\t' read -r SAMPLE CONDITION REP BAM; do
    SAMPLE_LIST+=("${SAMPLE}")
done < <(tail -n +2 "${SAMPLES}")

# --- Step 1: ATACorrect - parallel ------------------------------------------
echo "=== 1. ATACorrect - parallel ==="

PIDS=()
for SAMPLE in "${SAMPLE_LIST[@]}"; do
    OUT_BW="${TOBIAS_DIR}/${SAMPLE}.shifted_corrected.bw"
    if [[ -f "${OUT_BW}" ]]; then
        echo "  ${SAMPLE}: exists - skipping."
        continue
    fi

    echo "  ${SAMPLE}: starting..."
    TOBIAS ATACorrect \
      --bam "${SHIFT_DIR}/${SAMPLE}.shifted.bam" \
      --genome "${GENOME_FA}" \
      --peaks "${CONSENSUS_BED}" \
      --blacklist "${BLACKLIST}" \
      --outdir "${TOBIAS_DIR}" \
      --cores 1 \
      > "${TOBIAS_DIR}/${SAMPLE}_atacorrect.log" 2>&1 &

    PIDS+=($!)
    echo "  ${SAMPLE}: PID=${PIDS[-1]}"
done

echo "  Waiting for completion..."
for PID in "${PIDS[@]}"; do
    wait "${PID}"
    echo "  PID=${PID} done"
done
echo "  ATACorrect done!"

echo ""
echo "  Output files:"
for SAMPLE in "${SAMPLE_LIST[@]}"; do
    OUT_BW="${TOBIAS_DIR}/${SAMPLE}.shifted_corrected.bw"
    [[ -f "${OUT_BW}" ]] && echo "    OK: ${SAMPLE}" || echo "    FAIL: ${SAMPLE}"
done

# --- Step 2: ScoreBigwig - parallel -----------------------------------------
echo ""
echo "=== 2. ScoreBigwig - parallel ==="

PIDS=()
for SAMPLE in "${SAMPLE_LIST[@]}"; do
    CORRECTED_BW="${TOBIAS_DIR}/${SAMPLE}.shifted_corrected.bw"
    SCORE_BW="${TOBIAS_DIR}/${SAMPLE}_footprint.bw"

    if [[ -f "${SCORE_BW}" ]]; then
        echo "  ${SAMPLE}: exists - skipping."
        continue
    fi

    if [[ ! -f "${CORRECTED_BW}" ]]; then
        echo "  ${SAMPLE}: corrected BW not found - skipping."
        continue
    fi

    echo "  ${SAMPLE}: starting..."
    TOBIAS ScoreBigwig \
      --signal "${CORRECTED_BW}" \
      --regions "${CONSENSUS_BED}" \
      --output "${SCORE_BW}" \
      --cores 1 \
      > "${TOBIAS_DIR}/${SAMPLE}_score.log" 2>&1 &

    PIDS+=($!)
    echo "  ${SAMPLE}: PID=${PIDS[-1]}"
done

echo "  Waiting for completion..."
for PID in "${PIDS[@]}"; do
    wait "${PID}"
    echo "  PID=${PID} done"
done
echo "  ScoreBigwig done!"

# --- Step 3: BINDetect - per condition --------------------------------------
echo ""
echo "=== 3. BINDetect ==="

DMSO_BW1="${TOBIAS_DIR}/DMSO_Rep1_footprint.bw"
DMSO_BW2="${TOBIAS_DIR}/DMSO_Rep2_footprint.bw"

PIDS=()
for COND in PY60 PRBJN PRBJNsorted DUX4sorted; do
    BINDETECT_OUT="${TOBIAS_DIR}/BINDetect_${COND}_vs_DMSO"

    if [[ -f "${BINDETECT_OUT}/bindetect_results.txt" ]]; then
        echo "  ${COND}: exists - skipping."
        continue
    fi

    COND_BWS=()
    for SAMPLE in "${SAMPLE_LIST[@]}"; do
        COND_NAME=$(grep "^${SAMPLE}" "${SAMPLES}" | awk '{print $2}')
        [[ "${COND_NAME}" == "${COND}" ]] && \
            COND_BWS+=("${TOBIAS_DIR}/${SAMPLE}_footprint.bw")
    done

    mkdir -p "${BINDETECT_OUT}"
    echo "  ${COND}: starting (${#COND_BWS[@]} samples)..."

    TOBIAS BINDetect \
      --motifs "${JASPAR_MEME}" \
      --signals "${DMSO_BW1}" "${DMSO_BW2}" "${COND_BWS[@]}" \
      --genome "${GENOME_FA}" \
      --peaks "${CONSENSUS_BED}" \
      --outdir "${BINDETECT_OUT}" \
      --cond_names "DMSO_Rep1" "DMSO_Rep2" "${COND}_Rep1" "${COND}_Rep2" \
      --cores 1 \
      > "${TOBIAS_DIR}/BINDetect_${COND}.log" 2>&1 &

    PIDS+=($!)
    echo "  ${COND}: PID=${PIDS[-1]}"
done

echo "  Waiting for completion..."
for PID in "${PIDS[@]}"; do
    wait "${PID}"
    echo "  PID=${PID} done"
done

echo ""
echo "Done."
echo "  BINDetect: ${TOBIAS_DIR}/BINDetect_*/bindetect_results.txt"
