#!/bin/bash
#SBATCH --job-name=t5-2
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=t5-2_dexseq_summary_%j.out
#SBATCH --error=t5-2_dexseq_summary_%j.err

# ----------------------------------------
# Summarize DEXSeq count ouputs.  
# 
# Usage:
#   sbatch t5-3_dexseq_summary.sh
# ----------------------------------------

set -euo pipefail

BASE_DIR="${SLURM_SUBMIT_DIR:-$PWD}"

DEXSEQ_DIR="${BASE_DIR}/data/5_DEXSeq"
ANNOT_DIR="${DEXSEQ_DIR}/annotation"
COUNT_DIR="${DEXSEQ_DIR}/counts"
QC_DIR="${DEXSEQ_DIR}/count_qc"

SAMPLE_LIST="${BASE_DIR}/data/sample_list.txt"

SUMMARY="${QC_DIR}/dexseq_count_summary.tsv"

mkdir -p "${QC_DIR}"

# script info 
SCRIPT="t5-2_dexseq_summary.sh"
echo "******************************"
echo "Script name: ${SCRIPT}"
echo "SLURM job ID: ${SLURM_JOB_ID:-local}"
echo "Start time: $(date)"
echo "******************************"
echo ""

# summarize count outputs for each sample
echo -e "sample\tcount_file_exists\tn_lines\tn_assigned_features\tsum_counts" > "${SUMMARY}"

while read -r SAMPLE; do

    COUNT_FILE="${COUNT_DIR}/${SAMPLE}.dexseq_counts.txt"

    if [[ ! -f "${COUNT_FILE}" ]]; then
        echo -e "${SAMPLE}\tNO\t0\t0\t0" >> "${SUMMARY}"
        continue
    fi

    N_LINES=$(wc -l < "${COUNT_FILE}")

    N_ASSIGNED=$(awk '$1 !~ /^_/ {n++} END{print n+0}' "${COUNT_FILE}")
    SUM_COUNTS=$(awk '$1 !~ /^_/ {s+=$2} END{print s+0}' "${COUNT_FILE}")

    echo -e "${SAMPLE}\tYES\t${N_LINES}\t${N_ASSIGNED}\t${SUM_COUNTS}" >> "${SUMMARY}"

done < "${SAMPLE_LIST}"

echo "DEXSeq count summary saved to: ${SUMMARY}"
echo "End time: $(date)"
echo ""

