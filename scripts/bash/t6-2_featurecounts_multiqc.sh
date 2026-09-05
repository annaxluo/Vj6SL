#!/bin/bash
#SBATCH --job-name=t6-2
#SBATCH --output=t6-2_featurecounts_multiqc_%j.out
#SBATCH --error=t6-2_featurecounts_multiqc_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

# ----------------------------------------
# Summarize featureCounts outputs using MultiQC
# 
# Usage:
#   sbatch t6-2_featureCounts_multiqc.sh
# ----------------------------------------

set -euo pipefail
module load conda/3-23.3.1
conda activate rnaseq

# paths 
BASE_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
FC_DIR="${BASE_DIR}/data/6_featureCounts"
MULTIQC_DIR="${FC_DIR}/multiqc"
MULTIQC_FILENAME="featureCounts_multiqc.html"

LOG_DIR="${FC_DIR}/logs"
MULTIQC_STDOUT="${LOG_DIR}/multiqc_featureCounts.log"
MULTIQC_STDERR="${LOG_DIR}/multiqc_featureCounts.err"

mkdir -p "${MULTIQC_DIR}" "${LOG_DIR}"

# script info 
SCRIPT="t6-2_featureCounts_multiqc.sh"
echo "******************************"
echo "Script name: ${SCRIPT}"
echo "SLURM job ID: ${SLURM_JOB_ID:-local}"
echo "Start time: $(date)"
echo "******************************"
echo ""

# Run MultiQC
multiqc "${FC_DIR}" \
    -o "${MULTIQC_DIR}" \
    -n "${MULTIQC_FILENAME}" \
	--force \
    > "${MULTIQC_STDOUT}" \
    2> "${MULTIQC_STDERR}"

echo "MultiQC completed."
echo "Report saved to: ${MULTIQC_DIR}/${MULTIQC_FILENAME}"
echo "End time: $(date)"
