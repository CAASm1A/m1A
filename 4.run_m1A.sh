#!/bin/bash
# ============================================================
# run_all.sh — Run full m1A-seq downstream analysis pipeline
#
# Each step is a self-contained R script.
# Run all steps:   bash run_all.sh
# Run single step: bash run_all.sh 3   (runs only step 03)
# ============================================================

set -euo pipefail

RSCRIPT="${RSCRIPT:-Rscript}"   # Override with: RSCRIPT=/path/to/Rscript bash run_all.sh

STEPS=(
    "01_prepare_data.R"
    "02_annotate_peaks.R"
    "03_genomic_distribution.R"
    "04_motif_analysis.R"
    "05_sequence_features.R"
    "06_site_statistics.R"
    "07_go_enrichment.R"
    "08_pca.R"
)

run_step() {
    local script="$1"
    echo "============================================"
    echo "Running: ${script}"
    echo "$(date '+%Y-%m-%d %H:%M:%S')"
    echo "============================================"
    "${RSCRIPT}" "${script}"
    echo "[Done] ${script}"
    echo ""
}

if [[ $# -eq 1 ]]; then
    # Run a single step by number (e.g., bash run_all.sh 3)
    idx=$(( $1 - 1 ))
    run_step "${STEPS[$idx]}"
else
    # Run all steps sequentially
    for script in "${STEPS[@]}"; do
        run_step "${script}"
    done
    echo "All steps complete."
fi