#!/bin/bash
# ============================================================
# m1A-seq Analysis Pipeline — Step 1: QC & Alignment
# Organism: Arabidopsis thaliana (TAIR10)
# Usage:
#   bash m1A_pipeline.sh 1   # Step 1: Trim Galore QC
#   bash m1A_pipeline.sh 2   # Step 2: STAR alignment & quantification
#   bash m1A_pipeline.sh 3   # Step 3: deepTools visualization
# ============================================================

set -euo pipefail

# ─── User-configurable paths ────────────────────────────────
DIR="${PWD}"

# STAR genome index directories
GENOME_REF="/path/to/TAIR10/star"          # Main genome index
RIBO_REF="/path/to/TAIR10/ribo"            # Ribosomal RNA index (for depletion check)

# Annotation
GTF="/path/to/TAIR10/genes/genes.gtf"

# Path to STAR binary (e.g., from CellRanger or standalone install)
STAR_BIN="/path/to/STAR/bin"

# Conda environment containing trim_galore, samtools, deepTools, subread
CONDA_BIN="/path/to/conda/envs/tools/bin"

# ─── Samples ────────────────────────────────────────────────
# Edit this list to match your sample prefixes (must match raw FASTQ filenames)
# Expected filename pattern: ${sample}_R1.fastq.gz / ${sample}_R2.fastq.gz
SAMPLES="CIP1 H_IPDe1 I_IP2 J_IPDe1 K_IP1 K_IP2"

# ─── Parameters ─────────────────────────────────────────────
THREADS=20
SORT_THREADS=30
SORT_MEM="1G"
MAPQ=10          # Minimum mapping quality for filtering
BIN_SIZE=10      # deepTools bigwig bin size (bp)
BODY_LEN=5000    # Scale-regions body length (bp)
FLANK=3000       # Upstream/downstream flanking region (bp)

# ============================================================
# Step 1 — Adapter trimming & quality control (Trim Galore)
# ============================================================
if [[ "${1:-}" -eq 1 ]]; then
    echo "[Step 1] Running Trim Galore on: ${SAMPLES}"
    cd "${DIR}/rawdata"

    for sample in ${SAMPLES}; do
        echo "  Trimming: ${sample}"
        "${CONDA_BIN}/trim_galore" \
            -j "${THREADS}" \
            --paired \
            --length 30 \
            -q 20 \
            ${sample}_R1.fastq.gz ${sample}_R2.fastq.gz
    done

    echo "[Step 1] Done."
fi

# ============================================================
# Step 2 — STAR alignment, BAM processing & featureCounts
# ============================================================
if [[ "${1:-}" -eq 2 ]]; then
    echo "[Step 2] Aligning samples: ${SAMPLES}"

    # ── 2a. STAR alignment ──────────────────────────────────
    for sample in ${SAMPLES}; do
        FQ1="${DIR}/rawdata/${sample}_R1_val_1.fq.gz"
        FQ2="${DIR}/rawdata/${sample}_R2_val_2.fq.gz"

        if [[ ! -f "${FQ1}" || ! -f "${FQ2}" ]]; then
            echo "  [WARN] Trimmed FASTQs not found for ${sample}, skipping."
            continue
        fi

        echo "  Aligning: ${sample}"
        mkdir -p "${DIR}/matching/${sample}"

        "${STAR_BIN}/STAR" \
            --genomeDir "${GENOME_REF}" \
            --genomeLoad LoadAndKeep \
            --runThreadN "${THREADS}" \
            --readFilesIn "${FQ1}" "${FQ2}" \
            --readFilesCommand zcat \
            --outFileNamePrefix "${DIR}/matching/${sample}/${sample}_" \
            --outFilterMultimapNmax 20 \
            --outSAMtype BAM Unsorted \
            --outSAMstrandField intronMotif
    done

    # ── 2b. BAM sorting, filtering & bigwig generation ──────
    for sample in ${SAMPLES}; do
        OUT_DIR="${DIR}/matching/${sample}"
        RAW_BAM="${OUT_DIR}/${sample}_Aligned.out.bam"
        SORT_BAM="${OUT_DIR}/${sample}.sort.bam"
        UNIQ_BAM="${OUT_DIR}/${sample}.uniq.bam"
        BW="${OUT_DIR}/${sample}.bw"

        if [[ ! -f "${RAW_BAM}" ]]; then
            echo "  [WARN] BAM not found for ${sample}: ${RAW_BAM}"
            continue
        fi

        echo "  Processing BAM: ${sample}"

        # Sort
        samtools sort -@ "${SORT_THREADS}" -m "${SORT_MEM}" \
            "${RAW_BAM}" > "${SORT_BAM}"

        # Filter: unique mappers, properly paired (flags -q MAPQ -f 3)
        samtools view -b -q "${MAPQ}" -f 3 \
            "${SORT_BAM}" > "${UNIQ_BAM}"

        samtools index "${UNIQ_BAM}"

        # Generate normalised bigwig (BPM = reads per million mapped reads per bin)
        bamCoverage \
            --bam "${UNIQ_BAM}" \
            -o "${BW}" \
            --binSize "${BIN_SIZE}" \
            -p "${THREADS}" \
            --normalizeUsing BPM \
            --extendReads
    done

    # ── 2c. Gene-level quantification (featureCounts) ───────
    echo "  Running featureCounts..."
    ALL_BAMS=$(ls "${DIR}"/matching/*/*uniq.bam 2>/dev/null)

    featureCounts \
        -T "${SORT_THREADS}" \
        -a "${GTF}" \
        -o "${DIR}/all.featureCount.res" \
        ${ALL_BAMS}

    # ── 2d. FPKM calculation ─────────────────────────────────
    if [[ -f "${DIR}/get_FPKM.r" ]]; then
        echo "  Computing FPKM..."
        Rscript "${DIR}/get_FPKM.r"
    else
        echo "  [WARN] get_FPKM.r not found, skipping FPKM step."
    fi

    echo "[Step 2] Done."
fi

# ============================================================
# Step 3 — deepTools: matrix & profile plots
# ============================================================
if [[ "${1:-}" -eq 3 ]]; then
    echo "[Step 3] Running deepTools visualisation"

    mkdir -p "${DIR}/deeptools"
    cd "${DIR}/deeptools"

    # Symlink bigwig files for convenience
    ln -sf "${DIR}/matching/"*/*.bw ./

    for sample in ${SAMPLES}; do
        BW="${sample}.bw"
        if [[ ! -f "${BW}" ]]; then
            echo "  [WARN] Bigwig not found for ${sample}, skipping."
            continue
        fi

        echo "  Plotting: ${sample}"

        computeMatrix scale-regions \
            -p "${THREADS}" \
            -R "${GTF}" \
            -S "${BW}" \
            -b "${FLANK}" \
            -a "${FLANK}" \
            --regionBodyLength "${BODY_LEN}" \
            --skipZeros \
            -o "${sample}.gz"

        plotProfile \
            -m "${sample}.gz" \
            -out "${sample}.profiles.pdf" \
            --plotHeight 10 \
            --plotWidth 12 \
            --plotFileFormat pdf
    done

    echo "[Step 3] Done."
fi