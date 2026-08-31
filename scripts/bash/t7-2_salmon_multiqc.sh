#!/bin/bash
#SBATCH --job-name=t7-2
#SBATCH --output=t7-2_salmon_multiqc_%j.out
#SBATCH --error=t7-2_salmon_multiqc_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

# ----------------------------------------
# MultiQC summary for Salmon outputs
# 
# Usage:
#   sbatch t7-3_salmon_multiqc.sh
# ----------------------------------------

set -euo pipefail
module load conda/3-23.3.1
conda activate rnaseq

# paths 
BASE_DIR="${SLURM_SUBMIT_DIR:-$PWD}"

SALMON_DIR="${BASE_DIR}/data/7_Salmon"
MULTIQC_DIR="${SALMON_DIR}/multiqc"
MULTIQC_FILENAME="salmon_multiqc.html"
LOG_DIR="${SALMON_DIR}/logs"

mkdir -p "${MULTIQC_DIR}" "${LOG_DIR}" 

# script info 
ARRAY_TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
THREADS="${SLURM_CPUS_PER_TASK:-1}"
SCRIPT="t7-2_salmon_multiqc.sh"
echo "******************************"
echo "Script name: ${SCRIPT}"
echo "SLURM job ID: ${SLURM_JOB_ID:-local}"
echo "Start time: $(date)"
echo "******************************"
echo ""

multiqc "${SALMON_DIR}" \
    -o "${MULTIQC_DIR}" \
    -n "${MULTIQC_FILENAME}" \
	--force \
    > "${LOG_DIR}/multiqc_salmon.log" \
    2> "${LOG_DIR}/multiqc_salmon.err"

echo
echo "MultiQC completed."
echo "MultiQC report saved to: ${MULTIQC_DIR}/${MULTIQC_FILENAME}"
echo "End time: $(date)"
