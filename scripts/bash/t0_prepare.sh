#!/bin/bash
#SBATCH --job-name=t0
#SBATCH --output=t0_prepare_samples_%j.out
#SBATCH --error=t0_prepare_samples_%j.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00

# ---------------------------------------
# Prepare files for bulk RNA sequencing data processing workflow
# 1. Obtain a list of samples from raw data dir
# 2. Verify that each paired-end sample has both Read 1 and Read 2
#
# Usage:
#   sbatch t0_prepare.sh
# ---------------------------------------

set -euo pipefail

# paths 
BASE_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
INPUT_DIR="${BASE_DIR}/data/0_RawData"
OUTPUT_DIR="${BASE_DIR}/data"
SAMPLE_LIST="${OUTPUT_DIR}/sample_list.txt"

mkdir -p "${OUTPUT_DIR}"

# script info 
SCRIPT="t0_prepare.sh"
echo "******************************"
echo "Script name: ${SCRIPT}"
echo "SLURM job ID: ${SLURM_JOB_ID:-local}"
echo "Start date and time: $(date)"
echo "******************************"
echo ""

# 1. get sample list
echo "Creating sample list and checking fastq files..."
echo "Start time: $(date)"

ls ${INPUT_DIR}/*_1.fq.gz 2>/dev/null | \
    xargs -n1 basename | \
    sed 's/_1.fq.gz//' | \
    sort -V > "${SAMPLE_LIST}"

# 2. verify Read 1 and Read 2 file for each sample
echo "Verifying paired-end reads..."
echo "Start time: $(date)"
N_SAMPLES=$(wc -l < "${SAMPLE_LIST}")

while read -r SAMPLE; do
    R1="${INPUT_DIR}/${SAMPLE}_1.fq.gz"
    R2="${INPUT_DIR}/${SAMPLE}_2.fq.gz"

    if [[ -f "${R1}" && -f "${R2}" ]]; then
        echo "OK: ${SAMPLE}"
    else
        echo "Error: Missing pair for sample: ${SAMPLE}" >&2

        if [[ ! -f "${R1}" ]]; then
            echo "Missing Read 1: ${R1}" >&2
        fi

        if [[ ! -f "${R2}" ]]; then
            echo "Missing Read 2: ${R2}" >&2
        fi
    fi
done < "${SAMPLE_LIST}"

echo ""
echo "Sample list saved to: ${SAMPLE_LIST}"
echo "Script end time: $(date)"

