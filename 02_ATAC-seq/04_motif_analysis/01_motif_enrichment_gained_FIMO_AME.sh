#!/usr/bin/env bash
# =========================================================
# motif_enrichment_gained.sh
# ---------------------------------------------------------
# TF motif enrichment in gained peaks using the MEME Suite.
# For each condition:
#   1. Extract gained-peak sequences (hg38.fa)
#   2. Use all consensus-peak sequences as background
#   3. AME: gained vs background enrichment
#   4. FIMO: specific motif locations
#
# Requirements:
#   - MEME Suite (fimo, ame, fasta-get-markov)
#   - bedtools
# =========================================================

set -euo pipefail

# MEME Suite paths (adjust to your installation)
export MEME_DATA_DIR="$HOME/meme/share/meme-5.5.5"
export PATH="$HOME/meme/bin:$HOME/meme/libexec/meme-5.5.5:$PATH"

BASE_DIR="path/to/ATAC_analysis"
PROC_DIR="${BASE_DIR}/processed_multiNorm"
DIFF_DIR="${PROC_DIR}/diff_atac"
MOTIF_DIR="${PROC_DIR}/motif_analysis"
GENOME_FA="path/to/hg38/hg38.fa"
CONSENSUS_BED="${PROC_DIR}/counts/consensus_peaks.bed"
THREADS=6

mkdir -p "${MOTIF_DIR}"

# JASPAR2024 motif database
JASPAR_MEME="${MOTIF_DIR}/JASPAR2024_vertebrates.meme"
if [[ ! -f "${JASPAR_MEME}" ]]; then
    echo "Downloading JASPAR2024..."
    wget -q \
      "https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt" \
      -O "${JASPAR_MEME}"
fi

# --- Background: all consensus-peak sequences -------------------------------
echo "=== Preparing background sequences ==="

BG_FA="${MOTIF_DIR}/consensus_peaks_background.fa"
BG_MARKOV="${MOTIF_DIR}/background.markov"

if [[ ! -f "${BG_FA}" ]]; then
    echo "  Extracting consensus-peak sequences..."
    bedtools getfasta \
      -fi "${GENOME_FA}" \
      -bed "${CONSENSUS_BED}" \
      -fo "${BG_FA}"
fi

if [[ ! -f "${BG_MARKOV}" ]]; then
    echo "  Building Markov background model..."
    fasta-get-markov -m 2 "${BG_FA}" > "${BG_MARKOV}"
fi

# --- Motif analysis per condition -------------------------------------------
for COND in PY60 PRBJN PRBJNsorted DUX4sorted; do

    echo ""
    echo "================================"
    echo " ${COND} vs DMSO - motif analysis"
    echo "================================"

    GAINED_BED="${DIFF_DIR}/${COND}_vs_DMSO_gained.bed"

    if [[ ! -s "${GAINED_BED}" ]]; then
        echo "  BED empty - skipping."
        continue
    fi

    N=$(wc -l < "${GAINED_BED}")
    echo "  ${N} gained peaks"

    COND_DIR="${MOTIF_DIR}/${COND}"
    mkdir -p "${COND_DIR}"

    # --- Gained-peak FASTA ---
    GAINED_FA="${COND_DIR}/gained_peaks.fa"
    if [[ ! -f "${GAINED_FA}" ]]; then
        echo "  Extracting sequences..."
        bedtools getfasta \
          -fi "${GENOME_FA}" \
          -bed "${GAINED_BED}" \
          -fo "${GAINED_FA}"
    fi

    # --- AME: gained vs background enrichment ---
    AME_OUT="${COND_DIR}/ame_results"

    if [[ -d "${AME_OUT}" ]]; then
        echo "  AME results exist - skipping."
    else
        echo "  Running AME enrichment analysis... (${N} peaks)"
        ame \
          --oc "${AME_OUT}" \
          --control "${BG_FA}" \
          --method fisher \
          --scoring avg \
          --motif-pseudo 0.1 \
          "${GAINED_FA}" \
          "${JASPAR_MEME}"
        echo "  AME done: ${AME_OUT}/ame.html"
    fi

    # --- FIMO: specific motif locations ---
    FIMO_OUT="${COND_DIR}/fimo_results"

    if [[ -d "${FIMO_OUT}" ]]; then
        echo "  FIMO results exist - skipping."
    else
        echo "  Running FIMO motif scan..."
        fimo \
          --oc "${FIMO_OUT}" \
          --thresh 1e-4 \
          --bfile "${BG_MARKOV}" \
          "${JASPAR_MEME}" \
          "${GAINED_FA}"
    fi

    # --- Top motifs summary ---
    echo ""
    echo "  Top 10 enriched motifs (AME):"
    if [[ -f "${AME_OUT}/ame.tsv" ]]; then
        awk '!/^#/ && NF>5' "${AME_OUT}/ame.tsv" \
          | sort -k7,7g \
          | head -10 \
          | awk '{printf "    %-20s p=%-12s\n", $2, $7}'
    fi

done

echo ""
echo "Done."
echo "  Results: ${MOTIF_DIR}/<condition>/ame_results/ame.html"
