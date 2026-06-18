#!/usr/bin/env bash
# =========================================================
# homer_motif_enrichment.sh
# ---------------------------------------------------------
# TF motif enrichment in gained/lost peaks using HOMER.
# For each condition, runs findMotifsGenome.pl:
#   - Known-motif enrichment (knownResults.html)
#   - De novo motif discovery (homerResults.html)
#   - GC-matched background generated automatically
#
# Requirements:
#   - HOMER (findMotifsGenome.pl)
#   - hg38 installed: perl configureHomer.pl -install hg38
# =========================================================

set -euo pipefail

export PATH="$PATH:path/to/homer/bin"

BASE_DIR="path/to/ATAC_analysis"
PROC_DIR="${BASE_DIR}/processed_multiNorm"
DIFF_DIR="${PROC_DIR}/diff_atac"
HOMER_DIR="${PROC_DIR}/homer_motifs"
THREADS=6

mkdir -p "${HOMER_DIR}"

# --- Motif analysis per condition -------------------------------------------
for COND in PY60 PRBJN PRBJNsorted DUX4sorted; do

    echo ""
    echo "================================"
    echo " ${COND} vs DMSO - HOMER motif analysis"
    echo "================================"

    for DIRECTION in gained lost; do

        BED="${DIFF_DIR}/${COND}_vs_DMSO_${DIRECTION}.bed"

        if [[ ! -s "${BED}" ]]; then
            echo "  ${COND}_vs_DMSO_${DIRECTION}: empty BED - skipping."
            continue
        fi

        N=$(wc -l < "${BED}")
        OUT="${HOMER_DIR}/${COND}_vs_DMSO_${DIRECTION}"

        if [[ -f "${OUT}/knownResults.html" ]]; then
            echo "  ${COND}_${DIRECTION}: HOMER results exist - skipping."
            continue
        fi

        echo "  ${COND}_${DIRECTION}: ${N} peaks"
        mkdir -p "${OUT}"

        findMotifsGenome.pl \
          "${BED}" \
          hg38 \
          "${OUT}" \
          -size given \
          -mask \
          -p ${THREADS} \
          -preparse

        echo "  Done: ${OUT}/knownResults.html"

    done
done

# --- Top motifs summary -----------------------------------------------------
echo ""
echo "================================"
echo " Top 10 enriched known motifs"
echo "================================"

for COND in PY60 PRBJN PRBJNsorted DUX4sorted; do
    TSV="${HOMER_DIR}/${COND}_vs_DMSO_gained/knownResults.txt"
    if [[ -f "${TSV}" ]]; then
        echo ""
        echo "  ${COND} gained:"
        # Columns: Name, Consensus, P-value, Log P-value, q-value, ...
        awk -F'\t' 'NR>1 {printf "    %-30s p=%s\n", $1, $3}' "${TSV}" | head -10
    fi
done

echo ""
echo "Done."
echo "  Results : ${HOMER_DIR}/<condition>/knownResults.html"
echo "  De novo : ${HOMER_DIR}/<condition>/homerResults.html"
