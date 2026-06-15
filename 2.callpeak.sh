#!/bin/bash
# ============================================================
# m1A-seq Pipeline — Step 2: Peak Calling (MACS2)
# Usage: bash run_callpeak.sh
#
# Calls peaks for each sample pair (Input vs IP) in parallel.
# Expects callpeak.R to be in the same directory as this script.
# ============================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BAM_DIR="/path/to/bams"         # Directory containing all .bam files
RSCRIPT="/path/to/Rscript"      # e.g., /usr/bin/Rscript or conda env Rscript

# Sample prefixes and replicate numbers to process
SAMPLES="J"
REPLICATES="1 2"

cd "${BAM_DIR}"

for rep in ${REPLICATES}; do
    for sample in ${SAMPLES}; do
        INPUT="${sample}Input${rep}"
        IP="${sample}IP${rep}"
        echo "[callpeak] Input=${INPUT}  IP=${IP}"
        "${RSCRIPT}" "${SCRIPT_DIR}/callpeak.R" "${INPUT}" "${IP}" &
    done
done

wait
echo "All peak calling jobs finished."