#!/bin/bash
# ============================================================
# m1A-seq Pipeline — Step 3: Mutation Signature Analysis
#
# For each sample, within MACS2-called peak regions:
#   1. mpileup      — pileup reads at peak sites (IP & IPDe)
#   2. pileup2acgt  — convert pileup to per-base counts (A/C/G/T)
#   3. awk filters  — compute mismatch ratios; flag IPDe editing sites
#   4. merge        — intersect IP vs IPDe, call high-confidence m1A sites
#
# Output per sample (in ANALYSIS_DIR/):
#   ${sample}.final.xls     — candidate m1A sites (tab-separated)
#   ${sample}.final.bed     — BED file of candidate sites
#
# Usage: bash step3_mutation.sh
# ============================================================

set -euo pipefail

# ─── User-configurable paths ────────────────────────────────
GENOME_FA="/path/to/TAIR10/fasta/genome.fa"   # Reference genome FASTA
BAM_DIR="/path/to/bams"                        # Directory with *.uniq.bam files
PEAK_DIR="/path/to/bams"                       # Directory with MACS2 output subdirs
ANALYSIS_DIR="/path/to/analysis"               # Output directory

# ─── Parameters ─────────────────────────────────────────────
SAMPLES="A B C D E F G H I J K L M N O P"
REPLICATES="1 2"

# Filtering thresholds
MIN_DEPTH=4          # Minimum read depth at a site (>MIN_DEPTH)
MIN_RATIO=0.8        # Minimum mismatch ratio in IP to call a site
MIN_DELTA=0.2        # Minimum ratio difference (IP minus IPDe) to call enrichment
IPDe_NOISE=0.5       # IPDe sites with A% or T% above this are flagged as background

# ─── Main loop ───────────────────────────────────────────────
mkdir -p "${ANALYSIS_DIR}"
cd "${ANALYSIS_DIR}"

for rep in ${REPLICATES}; do
    for sample in ${SAMPLES}; do

        IP="${sample}IP${rep}"
        IPDe="${sample}IPDe${rep}"
        PEAK_BED="${PEAK_DIR}/${IP}/${IP}_peaks.narrowPeak"

        IP_BAM="${BAM_DIR}/${IP}.uniq.bam"
        IPDe_BAM="${BAM_DIR}/${IPDe}.uniq.bam"

        # Validate inputs
        for f in "${PEAK_BED}" "${IP_BAM}" "${IPDe_BAM}"; do
            if [[ ! -f "${f}" ]]; then
                echo "[WARN] Missing file, skipping ${IP}: ${f}"
                continue 2
            fi
        done

        (
            echo "[${IP}] Step 3a: mpileup"
            samtools mpileup -l "${PEAK_BED}" -f "${GENOME_FA}" "${IP_BAM}"    > "${IP}.pileup"   &
            samtools mpileup -l "${PEAK_BED}" -f "${GENOME_FA}" "${IPDe_BAM}"  > "${IPDe}.pileup" &
            wait

            echo "[${IP}] Step 3b: pileup2acgt"
            sequenza-utils pileup2acgt -p "${IP}.pileup"   > "${IP}.acgt"   &
            sequenza-utils pileup2acgt -p "${IPDe}.pileup" > "${IPDe}.acgt" &
            wait

            echo "[${IP}] Step 3c: compute per-base mismatch ratios"
            # Columns in .acgt: chr, pos, ref, depth, A, C, G, T, ...
            # Output: chr, pos, depth, A_ratio, T_ratio, ref_base
            awk 'NR>1 {
                OFS="\t"
                depth = $5+$6+$7+$8
                print $1, $2, depth, $5/(depth+1), $8/(depth+1), $3
            }' "${IP}.acgt" > "${IP}.t" &

            # IPDe: same ratios, keep only sites with high A or T mismatch (background noise)
            awk 'NR>1 {
                OFS="\t"
                depth = $5+$6+$7+$8
                print $1, $2, depth, $5/(depth+1), $8/(depth+1)
            }' "${IPDe}.acgt" \
            | awk -v thr="${IPDe_NOISE}" '$4 > thr || $5 > thr' \
            > "${IPDe}.t" &
            wait

            echo "[${IP}] Step 3d: merge IP & IPDe, filter candidate m1A sites"
            # Join on chr+pos: IPDe columns first, then IP columns
            awk -v ipde="${IPDe}.t" '
                BEGIN { while ((getline < ipde) > 0) lst[$1"\t"$2] = $0 }
                { key = $1"\t"$2; if (key in lst) print lst[key] "\t" $0 }
            ' "${IP}.t" > "${IP}_De_IP.res"

            # Filter: sufficient depth in both; IP mismatch ratio enriched vs IPDe
            # Columns after merge: IPDe(1-5) | IP(1-6)
            # $3=IPDe depth, $8=IP depth, $4=IPDe A_ratio, $9=IP A_ratio
            #                              $5=IPDe T_ratio, $10=IP T_ratio, $11=IP ref
            awk -v md="${MIN_DEPTH}" -v mr="${MIN_RATIO}" -v delta="${MIN_DELTA}" '
                $3 > md && $8 > md &&
                (($4 - $9 > delta && $4 > mr) || ($5 - $10 > delta && $5 > mr))
            ' "${IP}_De_IP.res" > "${IP}_De_IP.final.xls"

            # BED: chr, start, end (start+1), ref_base (col 11)
            awk 'OFS="\t" { print $1, $2, $2+1, $11 }' \
                "${IP}_De_IP.final.xls" > "${IP}_De_IP.final.bed"

            echo "[${IP}] Done. $(wc -l < "${IP}_De_IP.final.xls") candidate sites."
        ) &

    done
done

wait
echo "All samples processed."