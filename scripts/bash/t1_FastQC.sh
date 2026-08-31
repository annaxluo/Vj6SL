#!/bin/bash
#SBATCH --job-name=t1
#SBATCH --output=t1_fastqc_%j.out
#SBATCH --error=t1_fastqc_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00

# -------------------------------------
# Run pre-alignment QC
# 1. Perform FastQC for each raw fastq file
# 2. Summarize with MultiQC
#
# Usage:
#   sbatch t1_FastQC.sh
# -------------------------------------

set -euo pipefail

# modules
module load conda/3-23.3.1
conda activate seq3 
module load fastqc/0.12.1

# paths 
BASE_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
INPUT_DIR="${BASE_DIR}/data/0_RawData"
OUT_DIR="${BASE_DIR}/data/1_FastQC"
MULTIQC_DIR="${OUT_DIR}/MultiQC"

# Make output directory
mkdir -p "${OUT_DIR}"
mkdir -p "${MULTIQC_DIR}"

# script info 
SCRIPT="t1_FastQC.sh"
echo "******************************"
echo "Script name: ${SCRIPT}"
echo "SLURM job ID: ${SLURM_JOB_ID:-local}"
echo "Start date and time: $(date)"
echo "******************************"
echo ""

# 1. run FastQC
echo "Running FastQC..."
echo "Start time: $(date)"

fastqc \
    --threads "${SLURM_CPUS_PER_TASK:-1}" \
    --outdir "${OUT_DIR}" \
    "${INPUT_DIR}"/*.fq.gz

echo "FastQC complete."
echo "End time: $(date)"
echo ""

# 2. run MultiQC
echo "Running MultiQC..."

multiqc \
    "${OUT_DIR}" \
    --outdir "${MULTIQC_DIR}" \
    --filename "multiqc_fastqc_raw.html" \
    --title "Raw FASTQ FastQC Summary"

echo "FastQC + MultiQC complete."
echo "Script end time: $(date)"