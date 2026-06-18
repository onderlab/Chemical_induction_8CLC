#!/usr/bin/env bash
# =========================================================
# ATAC-seq main pipeline (BAM -> normalized signal tracks)
# =========================================================
# 12-step workflow:
#   1.  Filter + blacklist removal + ATAC shift
#   2.  Peak calling (MACS2)
#   3.  Consensus peaks (merge across samples)
#   4.  featureCounts on consensus peaks
#   5.  TSS BED files for all / ZGA / housekeeping genes
#   6.  R normalization (stable peaks + HK promoters)
#   7.  bigWig generation: CPM / RPGC / stablePeak
#   8.  computeMatrix + plotProfile + plotHeatmap (CPM/RPGC/stablePeak)
#   9.  Numeric comparison summaries across normalizations
#   10. hkTSS scale-factor computation (HK TSS read counts)
#   11. hkTSS bigWig generation
#   12. computeMatrix + plotProfile + plotHeatmap (hkTSS)
# ---------------------------------------------------------

set -euo pipefail

# Fix macOS open file limit for alignmentSieve
ulimit -n 10000

# ---- Paths (replace with your own) ----
BASE_DIR="path/to/ATAC_analysis"
ALIGNED_DIR="${BASE_DIR}/aligned"
PROC_DIR="${BASE_DIR}/processed_multiNorm"
SHIFT_DIR="${PROC_DIR}/shifted_bam"
PEAK_DIR="${PROC_DIR}/peaks"
COUNT_DIR="${PROC_DIR}/counts"
BW_DIR="${PROC_DIR}/bigwig"
MATRIX_DIR="${PROC_DIR}/matrices"
BED_DIR="${PROC_DIR}/beds"
QC_DIR="${PROC_DIR}/qc"
SCRIPT_DIR="${BASE_DIR}/scripts"

GENOME_DIR="path/to/hg38"
GTF="${GENOME_DIR}/hg38_genes.gtf"
CHROM_SIZES="${GENOME_DIR}/hg38.chrom.sizes"
BLACKLIST="${GENOME_DIR}/hg38-blacklist.v2.bed"

SAMPLES="${BASE_DIR}/samples.tsv"
ZGA_GENES="${BASE_DIR}/zga_genes.txt"
HK_GENES="${BASE_DIR}/housekeeping_genes.txt"

EFFECTIVE_GENOME_SIZE=2913022398
THREADS=6

mkdir -p "${PROC_DIR}" "${SHIFT_DIR}" "${PEAK_DIR}" "${COUNT_DIR}" "${BW_DIR}" \
         "${MATRIX_DIR}" "${BED_DIR}" "${QC_DIR}" "${SCRIPT_DIR}"

# Strip Windows \r from input files if present
sed -i '' 's/\r$//' "${SAMPLES}" 2>/dev/null || sed -i 's/\r$//' "${SAMPLES}"
sed -i '' 's/\r$//' "${ZGA_GENES}" 2>/dev/null || sed -i 's/\r$//' "${ZGA_GENES}"
sed -i '' 's/\r$//' "${HK_GENES}" 2>/dev/null || sed -i 's/\r$//' "${HK_GENES}"

# =========================================================
# STEP 1: Filter + blacklist removal + ATAC shift
# =========================================================
echo "=== 1. Filter + blacklist removal + ATAC shift ==="

STEP1_SKIP=true
tail -n +2 "${SAMPLES}" | while IFS=$'\t' read -r SAMPLE CONDITION REP BAM; do
    if [[ ! -f "${SHIFT_DIR}/${SAMPLE}.shifted.bam" || ! -f "${SHIFT_DIR}/${SAMPLE}.shifted.bam.bai" ]]; then
        exit 1
    fi
done && STEP1_SKIP=true || STEP1_SKIP=false

if ${STEP1_SKIP}; then
    echo "  All shifted BAMs exist - skipping Step 1."
else
    tail -n +2 "${SAMPLES}" | while IFS=$'\t' read -r SAMPLE CONDITION REP BAM; do

        if [[ -f "${SHIFT_DIR}/${SAMPLE}.shifted.bam" && -f "${SHIFT_DIR}/${SAMPLE}.shifted.bam.bai" ]]; then
            echo "  ${SAMPLE} already done - skipping."
            continue
        fi

        echo "  Processing ${SAMPLE}"

        TMP1="${SHIFT_DIR}/${SAMPLE}.tmp1.bam"
        TMP2="${SHIFT_DIR}/${SAMPLE}.tmp2.bam"
        FILT="${SHIFT_DIR}/${SAMPLE}.filtered.bam"
        SHIFT="${SHIFT_DIR}/${SAMPLE}.shifted.bam"

        # Quality filter: properly paired (-f 2), drop unmapped/secondary/dup (-F 1804), MAPQ >= 30
        # Remove mitochondrial reads
        samtools view -@ ${THREADS} -h -b -f 2 -F 1804 -q 30 "${BAM}" \
          | samtools view -@ ${THREADS} -h \
          | awk 'BEGIN{OFS="\t"} /^@/ {print; next} $3!="chrM" {print}' \
          | samtools sort -@ ${THREADS} -o "${TMP1}"

        samtools index "${TMP1}"

        # Remove ENCODE-blacklisted regions
        bedtools intersect -v -abam "${TMP1}" -b "${BLACKLIST}" > "${TMP2}"
        samtools sort -@ ${THREADS} -o "${FILT}" "${TMP2}"
        samtools index "${FILT}"

        # ATAC-seq Tn5 insertion-bias correction (+4 / -5 bp shift)
        alignmentSieve --bam "${FILT}" --ATACshift -p ${THREADS} \
          -o "${SHIFT_DIR}/${SAMPLE}.shifted.unsorted.bam"

        samtools sort -@ ${THREADS} -o "${SHIFT}" "${SHIFT_DIR}/${SAMPLE}.shifted.unsorted.bam"
        samtools index "${SHIFT}"
        rm -f "${SHIFT_DIR}/${SAMPLE}.shifted.unsorted.bam"

        samtools flagstat "${SHIFT}" > "${QC_DIR}/${SAMPLE}.flagstat.txt"

        rm -f "${TMP1}" "${TMP1}.bai" "${TMP2}"
    done
fi

# =========================================================
# STEP 2: Peak calling (MACS2)
# =========================================================
echo "=== 2. Peak calling ==="

STEP2_SKIP=true
tail -n +2 "${SAMPLES}" | while IFS=$'\t' read -r SAMPLE CONDITION REP BAM; do
    [[ ! -f "${PEAK_DIR}/${SAMPLE}_peaks.narrowPeak" ]] && exit 1
    true
done && STEP2_SKIP=true || STEP2_SKIP=false

if ${STEP2_SKIP}; then
    echo "  All narrowPeak files exist - skipping Step 2."
else
    tail -n +2 "${SAMPLES}" | while IFS=$'\t' read -r SAMPLE CONDITION REP BAM; do

        if [[ -f "${PEAK_DIR}/${SAMPLE}_peaks.narrowPeak" ]]; then
            echo "  ${SAMPLE} peaks exist - skipping."
            continue
        fi

        echo "  Calling peaks for ${SAMPLE}"
        SHIFT="${SHIFT_DIR}/${SAMPLE}.shifted.bam"

        macs2 callpeak \
          -t "${SHIFT}" \
          -f BAMPE \
          -g hs \
          -n "${SAMPLE}" \
          --outdir "${PEAK_DIR}" \
          --keep-dup all \
          -q 0.01
    done
fi

# =========================================================
# STEP 3: Consensus peaks (merge across all samples)
# =========================================================
echo "=== 3. Consensus peaks ==="

if [[ -f "${COUNT_DIR}/consensus_peaks.bed" && -f "${COUNT_DIR}/consensus_peaks.saf" ]]; then
    echo "  Consensus peaks exist - skipping Step 3."
else
    cat "${PEAK_DIR}"/*_peaks.narrowPeak \
      | awk 'BEGIN{OFS="\t"} $1 ~ /^chr([0-9]+|X|Y)$/ {print $1,$2,$3}' \
      | sort -k1,1 -k2,2n \
      | bedtools merge -i - \
      > "${COUNT_DIR}/consensus_peaks.bed"

    awk 'BEGIN{OFS="\t"; print "GeneID","Chr","Start","End","Strand"}
         {print "peak_"NR,$1,$2+1,$3,"."}' \
      "${COUNT_DIR}/consensus_peaks.bed" > "${COUNT_DIR}/consensus_peaks.saf"

    echo "  Consensus peaks: $(wc -l < "${COUNT_DIR}/consensus_peaks.bed") regions"
fi

# =========================================================
# STEP 4: Count consensus peaks (featureCounts)
# =========================================================
echo "=== 4. Count consensus peaks ==="

RERUN_STEP4=false
if [[ -f "${COUNT_DIR}/consensus_peak_counts.txt" ]]; then
    EXPECTED=$(tail -n +2 "${SAMPLES}" | wc -l | tr -d ' ')
    GOT=$(head -2 "${COUNT_DIR}/consensus_peak_counts.txt" | tail -1 | awk -F'\t' '{print NF - 6}')
    if [[ "${GOT}" -eq "${EXPECTED}" ]]; then
        echo "  Count file exists with ${GOT}/${EXPECTED} samples - skipping Step 4."
    else
        echo "  Count file has ${GOT} samples but expected ${EXPECTED} - re-running."
        RERUN_STEP4=true
    fi
else
    RERUN_STEP4=true
fi

if ${RERUN_STEP4}; then
    BAMS=()
    while IFS=$'\t' read -r SAMPLE CONDITION REP BAM; do
        [[ "${SAMPLE}" == "sample" ]] && continue
        BAMS+=("${SHIFT_DIR}/${SAMPLE}.shifted.bam")
    done < "${SAMPLES}"

    featureCounts \
      -a "${COUNT_DIR}/consensus_peaks.saf" \
      -F SAF \
      -T ${THREADS} \
      -p -B -C \
      -o "${COUNT_DIR}/consensus_peak_counts.txt" \
      "${BAMS[@]}"
fi

# =========================================================
# STEP 5: TSS BED files (all / ZGA / housekeeping)
# =========================================================
echo "=== 5. Build TSS BEDs for all / ZGA / housekeeping ==="

if [[ -f "${BED_DIR}/all_gene_tss_1bp.bed" && -f "${BED_DIR}/zga_tss_1bp.bed" && -f "${BED_DIR}/hk_tss_1bp.bed" ]]; then
    echo "  TSS BEDs exist - skipping Step 5."
else
    # Extract TSS coordinates from GTF (1-bp interval at gene start, strand-aware)
    gawk 'BEGIN{FS=OFS="\t"}
    $3=="gene"{
        gene_name=""; strand=$7;
        if (match($9, /gene_name "([^"]+)"/, a)) gene_name=a[1];
        if (gene_name!="") {
            if (strand=="+") {tss=$4-1} else {tss=$5-1}
            print $1,tss,tss+1,gene_name,".",strand
        }
    }' "${GTF}" | awk '$1 ~ /^chr([0-9]+|X|Y)$/' > "${BED_DIR}/all_gene_tss_1bp.bed"

    # Subset to ZGA and housekeeping gene lists
    awk 'NR==FNR {keep[$1]=1; next} ($4 in keep)' "${ZGA_GENES}" "${BED_DIR}/all_gene_tss_1bp.bed" \
      > "${BED_DIR}/zga_tss_1bp.bed"

    awk 'NR==FNR {keep[$1]=1; next} ($4 in keep)' "${HK_GENES}" "${BED_DIR}/all_gene_tss_1bp.bed" \
      > "${BED_DIR}/hk_tss_1bp.bed"

    echo "  all_gene_tss: $(wc -l < "${BED_DIR}/all_gene_tss_1bp.bed")"
    echo "  zga_tss:      $(wc -l < "${BED_DIR}/zga_tss_1bp.bed")"
    echo "  hk_tss:       $(wc -l < "${BED_DIR}/hk_tss_1bp.bed")"
fi

# =========================================================
# STEP 6: R normalization (stable peaks + HK promoters)
# =========================================================
echo "=== 6. R normalization (stable peaks + HK promoters) ==="

if [[ -f "${COUNT_DIR}/stablePeak_scale_factors.tsv" && -f "${QC_DIR}/PCA_stablePeakNorm.pdf" ]]; then
    echo "  Scale factors and PCA exist - skipping Step 6."
else
    Rscript "${SCRIPT_DIR}/atac_multinorm.R" \
      "${COUNT_DIR}/consensus_peak_counts.txt" \
      "${SAMPLES}" \
      "${BED_DIR}/hk_tss_1bp.bed" \
      "${COUNT_DIR}" \
      "${QC_DIR}"
fi

# =========================================================
# STEP 7: bigWig generation (CPM / RPGC / stablePeak)
# =========================================================
echo "=== 7. bigWig generation: CPM / RPGC / stablePeak ==="

# CPM + RPGC
STEP7_CPMRPGC_SKIP=true
tail -n +2 "${SAMPLES}" | while IFS=$'\t' read -r SAMPLE CONDITION REP BAM; do
    [[ ! -f "${BW_DIR}/${SAMPLE}.CPM.bw" || ! -f "${BW_DIR}/${SAMPLE}.RPGC.bw" ]] && exit 1
    true
done && STEP7_CPMRPGC_SKIP=true || STEP7_CPMRPGC_SKIP=false

if ${STEP7_CPMRPGC_SKIP}; then
    echo "  All CPM and RPGC bigWigs exist - skipping."
else
    tail -n +2 "${SAMPLES}" | while IFS=$'\t' read -r SAMPLE CONDITION REP BAM; do
        SHIFT="${SHIFT_DIR}/${SAMPLE}.shifted.bam"

        if [[ ! -f "${BW_DIR}/${SAMPLE}.CPM.bw" ]]; then
            echo "  CPM bigWig: ${SAMPLE}"
            bamCoverage -b "${SHIFT}" -o "${BW_DIR}/${SAMPLE}.CPM.bw" \
              --binSize 10 --normalizeUsing CPM --extendReads -p ${THREADS}
        fi

        if [[ ! -f "${BW_DIR}/${SAMPLE}.RPGC.bw" ]]; then
            echo "  RPGC bigWig: ${SAMPLE}"
            bamCoverage -b "${SHIFT}" -o "${BW_DIR}/${SAMPLE}.RPGC.bw" \
              --binSize 10 --normalizeUsing RPGC \
              --effectiveGenomeSize "${EFFECTIVE_GENOME_SIZE}" --extendReads -p ${THREADS}
        fi
    done
fi

# stablePeak (scaled by stable-peak size factors from Step 6)
STEP7_SP_SKIP=true
if [[ -f "${COUNT_DIR}/stablePeak_scale_factors.tsv" ]]; then
    tail -n +2 "${COUNT_DIR}/stablePeak_scale_factors.tsv" | while IFS=$'\t' read -r SAMPLE SIZEFACTOR SCALEFACTOR; do
        [[ ! -f "${BW_DIR}/${SAMPLE}.stablePeak.bw" ]] && exit 1
        true
    done && STEP7_SP_SKIP=true || STEP7_SP_SKIP=false
else
    STEP7_SP_SKIP=false
fi

if ${STEP7_SP_SKIP}; then
    echo "  All stablePeak bigWigs exist - skipping."
else
    tail -n +2 "${COUNT_DIR}/stablePeak_scale_factors.tsv" | while IFS=$'\t' read -r SAMPLE SIZEFACTOR SCALEFACTOR; do
        SHIFT="${SHIFT_DIR}/${SAMPLE}.shifted.bam"

        if [[ ! -f "${BW_DIR}/${SAMPLE}.stablePeak.bw" ]]; then
            echo "  stablePeak bigWig: ${SAMPLE}"
            bamCoverage -b "${SHIFT}" -o "${BW_DIR}/${SAMPLE}.stablePeak.bw" \
              --binSize 10 --normalizeUsing None \
              --scaleFactor "${SCALEFACTOR}" --extendReads -p ${THREADS}
        fi
    done
fi

# =========================================================
# STEP 8: computeMatrix + plotProfile + plotHeatmap
# =========================================================
echo "=== 8. computeMatrix for all TSS / HK TSS / ZGA TSS ==="

SAMPLE_LABELS=$(cut -f1 "${SAMPLES}" | tail -n +2 | tr '\n' ' ')

for MODE in CPM RPGC stablePeak; do

    # Check if ALL outputs for this MODE exist
    ALL_DONE=true
    for REG in allTSS hkTSS zgaTSS; do
        [[ ! -f "${MATRIX_DIR}/${REG}_${MODE}.gz" ]] && ALL_DONE=false
        [[ ! -f "${MATRIX_DIR}/${REG}_${MODE}_profile.pdf" ]] && ALL_DONE=false
        [[ ! -f "${MATRIX_DIR}/${REG}_${MODE}_profile.tsv" ]] && ALL_DONE=false
    done
    [[ ! -f "${MATRIX_DIR}/zgaTSS_${MODE}_heatmap.pdf" ]] && ALL_DONE=false

    if ${ALL_DONE}; then
        echo "  ${MODE}: all outputs exist - skipping."
        continue
    fi

    # Build bigWig file list
    if [[ "${MODE}" == "stablePeak" ]]; then
        BWFILES=()
        LABELS=()
        while IFS=$'\t' read -r S _SF _SC; do
            BWFILES+=("${BW_DIR}/${S}.${MODE}.bw")
            LABELS+=("${S}")
        done < <(tail -n +2 "${COUNT_DIR}/stablePeak_scale_factors.tsv")
        LABEL_STR="${LABELS[*]}"
    else
        BWFILES=()
        while IFS=$'\t' read -r SAMPLE CONDITION REP BAM; do
            [[ "${SAMPLE}" == "sample" ]] && continue
            BWFILES+=("${BW_DIR}/${SAMPLE}.${MODE}.bw")
        done < "${SAMPLES}"
        LABEL_STR="${SAMPLE_LABELS}"
    fi

    for REG in allTSS hkTSS zgaTSS; do
        case "${REG}" in
            allTSS) BED="${BED_DIR}/all_gene_tss_1bp.bed" ;;
            hkTSS)  BED="${BED_DIR}/hk_tss_1bp.bed" ;;
            zgaTSS) BED="${BED_DIR}/zga_tss_1bp.bed" ;;
        esac

        if [[ ! -f "${MATRIX_DIR}/${REG}_${MODE}.gz" ]]; then
            echo "  Computing matrix: ${REG}_${MODE}"
            computeMatrix reference-point \
              --referencePoint TSS -b 2000 -a 2000 \
              -R "${BED}" -S "${BWFILES[@]}" \
              --skipZeros --missingDataAsZero -p ${THREADS} \
              -o "${MATRIX_DIR}/${REG}_${MODE}.gz"
        fi

        if [[ ! -f "${MATRIX_DIR}/${REG}_${MODE}_profile.pdf" || ! -f "${MATRIX_DIR}/${REG}_${MODE}_profile.tsv" ]]; then
            echo "  Plotting profile: ${REG}_${MODE}"
            plotProfile \
              -m "${MATRIX_DIR}/${REG}_${MODE}.gz" \
              --perGroup \
              --samplesLabel ${LABEL_STR} \
              --refPointLabel "TSS" \
              --outFileNameData "${MATRIX_DIR}/${REG}_${MODE}_profile.tsv" \
              -out "${MATRIX_DIR}/${REG}_${MODE}_profile.pdf"
        fi
    done

    if [[ ! -f "${MATRIX_DIR}/zgaTSS_${MODE}_heatmap.pdf" ]]; then
        echo "  Plotting heatmap: zgaTSS_${MODE}"
        plotHeatmap -m "${MATRIX_DIR}/zgaTSS_${MODE}.gz" \
          --colorMap Reds \
          --whatToShow 'heatmap and colorbar' \
          --samplesLabel ${LABEL_STR} \
          --refPointLabel "TSS" \
          --xAxisLabel "Distance from TSS (bp)" \
          --zMin 0 \
          --sortUsing mean \
          --sortRegions descend \
          --missingDataColor '#f7f7f7' \
          --heatmapHeight 15 \
          --heatmapWidth 3 \
          -out "${MATRIX_DIR}/zgaTSS_${MODE}_heatmap.pdf"
    fi
done

# =========================================================
# STEP 9: Numeric summaries (compare normalizations)
# =========================================================
echo "=== 9. Extract numeric summaries for direct comparison ==="

if [[ -f "${QC_DIR}/normalization_comparison_summary.csv" && -f "${QC_DIR}/HK_TSS_CV_by_normalization.csv" ]]; then
    echo "  Comparison summaries exist - skipping Step 9."
else
    Rscript "${SCRIPT_DIR}/compare_region_sets.R" \
      "${MATRIX_DIR}" \
      "${SAMPLES}" \
      "${QC_DIR}"
fi

# =========================================================
# STEP 10: hkTSS scale-factor computation
# =========================================================
echo "=== 10. hkTSS scale factor computation ==="

if [[ -f "${COUNT_DIR}/hkTSS_scale_factors.tsv" ]]; then
    echo "  hkTSS scale factors exist - skipping."
else
    Rscript "${SCRIPT_DIR}/hkTSS_normalize.R" \
      "${BED_DIR}/hk_tss_1bp.bed" \
      "${SHIFT_DIR}" \
      "${SAMPLES}" \
      "${COUNT_DIR}" \
      "${QC_DIR}"
fi

# =========================================================
# STEP 11: hkTSS bigWig generation
# =========================================================
echo "=== 11. hkTSS bigWig generation ==="

STEP11_SKIP=true
tail -n +2 "${COUNT_DIR}/hkTSS_scale_factors.tsv" | while IFS=$'\t' read -r SAMPLE _SF SC; do
    [[ ! -f "${BW_DIR}/${SAMPLE}.hkTSS.bw" ]] && exit 1
    true
done && STEP11_SKIP=true || STEP11_SKIP=false

if ${STEP11_SKIP}; then
    echo "  All hkTSS bigWigs exist - skipping."
else
    tail -n +2 "${COUNT_DIR}/hkTSS_scale_factors.tsv" | while IFS=$'\t' read -r SAMPLE _SF SC; do
        SHIFT="${SHIFT_DIR}/${SAMPLE}.shifted.bam"
        OUT="${BW_DIR}/${SAMPLE}.hkTSS.bw"

        if [[ -f "${OUT}" ]]; then
            echo "  ${SAMPLE} hkTSS bigWig exists - skipping."
            continue
        fi

        echo "  hkTSS bigWig: ${SAMPLE} (scaleFactor=${SC})"
        bamCoverage -b "${SHIFT}" -o "${OUT}" \
          --binSize 10 --normalizeUsing None \
          --scaleFactor "${SC}" --extendReads -p ${THREADS}
    done
fi

# =========================================================
# STEP 12: computeMatrix + plotProfile + plotHeatmap (hkTSS)
# =========================================================
echo "=== 12. computeMatrix for hkTSS normalisation ==="

BWFILES_HK=()
LABELS_HK=()
while IFS=$'\t' read -r SAMPLE _SF _SC; do
    BWFILES_HK+=("${BW_DIR}/${SAMPLE}.hkTSS.bw")
    LABELS_HK+=("${SAMPLE}")
done < <(tail -n +2 "${COUNT_DIR}/hkTSS_scale_factors.tsv")
LABEL_STR_HK="${LABELS_HK[*]}"

ALL_DONE_HK=true
for REG in allTSS hkTSS zgaTSS; do
    [[ ! -f "${MATRIX_DIR}/${REG}_hkTSS.gz" ]]          && ALL_DONE_HK=false
    [[ ! -f "${MATRIX_DIR}/${REG}_hkTSS_profile.pdf" ]] && ALL_DONE_HK=false
    [[ ! -f "${MATRIX_DIR}/${REG}_hkTSS_profile.tsv" ]] && ALL_DONE_HK=false
done
[[ ! -f "${MATRIX_DIR}/zgaTSS_hkTSS_heatmap.pdf" ]] && ALL_DONE_HK=false

if ${ALL_DONE_HK}; then
    echo "  hkTSS: all outputs exist - skipping."
else
    for REG in allTSS hkTSS zgaTSS; do
        case "${REG}" in
            allTSS) BED="${BED_DIR}/all_gene_tss_1bp.bed" ;;
            hkTSS)  BED="${BED_DIR}/hk_tss_1bp.bed" ;;
            zgaTSS) BED="${BED_DIR}/zga_tss_1bp.bed" ;;
        esac

        if [[ ! -f "${MATRIX_DIR}/${REG}_hkTSS.gz" ]]; then
            echo "  Computing matrix: ${REG}_hkTSS"
            computeMatrix reference-point \
              --referencePoint TSS -b 2000 -a 2000 \
              -R "${BED}" -S "${BWFILES_HK[@]}" \
              --skipZeros --missingDataAsZero -p ${THREADS} \
              -o "${MATRIX_DIR}/${REG}_hkTSS.gz"
        fi

        if [[ ! -f "${MATRIX_DIR}/${REG}_hkTSS_profile.pdf" || ! -f "${MATRIX_DIR}/${REG}_hkTSS_profile.tsv" ]]; then
            echo "  Plotting profile: ${REG}_hkTSS"
            plotProfile \
              -m "${MATRIX_DIR}/${REG}_hkTSS.gz" \
              --perGroup \
              --samplesLabel ${LABEL_STR_HK} \
              --refPointLabel "TSS" \
              --outFileNameData "${MATRIX_DIR}/${REG}_hkTSS_profile.tsv" \
              -out "${MATRIX_DIR}/${REG}_hkTSS_profile.pdf"
        fi
    done

    if [[ ! -f "${MATRIX_DIR}/zgaTSS_hkTSS_heatmap.pdf" ]]; then
        echo "  Plotting heatmap: zgaTSS_hkTSS"
        plotHeatmap -m "${MATRIX_DIR}/zgaTSS_hkTSS.gz" \
          --colorMap Reds \
          --whatToShow 'heatmap and colorbar' \
          --samplesLabel ${LABEL_STR_HK} \
          --refPointLabel "TSS" \
          --xAxisLabel "Distance from TSS (bp)" \
          --zMin 0 \
          --sortUsing mean \
          --sortRegions descend \
          --missingDataColor '#f7f7f7' \
          --heatmapHeight 15 \
          --heatmapWidth 3 \
          -out "${MATRIX_DIR}/zgaTSS_hkTSS_heatmap.pdf"
    fi
fi

echo "Done."
