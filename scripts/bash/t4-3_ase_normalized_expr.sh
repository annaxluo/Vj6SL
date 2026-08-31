#!/bin/bash
#SBATCH --job-name=t4-3
#SBATCH --output=t4-3_ase_normalized_expr_%j.out
#SBATCH --error=t4-3_ase_normalized_expr_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

# ----------------------------------------
# Normalize REF and ALT counts to total library size. 
# 1. Count mapped reads per sample. 
# 2. Compute normalized REF and ALT counts. 
# 
# Usage:
#   sbatch t4-4_ase_normalized_expr.sh
#   sbatch t4-4_ase_normalized_expr.sh 10
# ----------------------------------------

set -eo pipefail

# input validity 
MIN_DEPTH="${1:-10}"

if ! [[ "${MIN_DEPTH}" =~ ^[0-9]+$ ]]; then
    echo "Error: MIN_DEPTH must be a positive integer."
    exit 1
fi

set +u
module load samtools/1.18
module load conda/3-23.3.1
conda activate rnaseq
set -u

# paths
BASE_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
ALIGN_DIR="${BASE_DIR}/data/3_STAR_Alignment"
FILTERED_BAM_DIR="${ALIGN_DIR}/filtered_bams"

ASE_DIR="${BASE_DIR}/data/4_ASE_Analysis"
ASE_COUNT_DIR="${ASE_DIR}/allele_counts_Q20"

SAMPLE_LIST="${BASE_DIR}/data/sample_list.txt"
SAMPLESHEET="${BASE_DIR}/data/samplesheet.csv"

OUT_DIR="${ASE_DIR}/normalized_expression"
MAPPED_READS="${OUT_DIR}/mapped_reads.tsv"
OUT_TSV="${OUT_DIR}/ASE_expression.normalized.tsv"

PYTHON_SCRIPT="${BASE_DIR}/python/t4_normalize_expression.py"

mkdir -p "${OUT_DIR}"

# script info 
SCRIPT="t4-3_ase_normalized_expr.sh"
echo "******************************"
echo "Script name: ${SCRIPT}"
echo "SLURM job ID: ${SLURM_JOB_ID:-local}"
echo "Start date and time: $(date)"
echo "******************************"
echo ""

# sample list 
if [[ ! -f "${SAMPLE_LIST}" ]]; then
    echo "Error: sample list not found: ${SAMPLE_LIST}" 
    exit 1
fi

mapfile -t SAMPLES < "${SAMPLE_LIST}"
N_SAMPLES=$(wc -l < "${SAMPLE_LIST}")

echo "Sample list: ${SAMPLE_LIST}"
echo "Total samples in list: ${N_SAMPLES}"
echo ""


# 1. Computing mapped reads from read-group BAM files
echo "1. Computing mapped reads...."
echo "Start time: $(date)"

echo -e "sample_id\tmapped_reads" > "${MAPPED_READS}"

while read -r SAMPLE; do
    SAMPLE=$(echo "${SAMPLE}" | tr -d '\r')

    BAM="${FILTERED_BAM_DIR}/${SAMPLE}.wasp_filtered.rg.bam"

    MAPPED=$(samtools idxstats "${BAM}" | awk '{sum += $3} END {print sum + 0}')

    echo -e "${SAMPLE}\t${MAPPED}" >> "${MAPPED_READS}"

done < "${SAMPLE_LIST}"

echo "Mapped reads written to: ${MAPPED_READS}"
echo "End time: $(date)"
echo ""

# 2. Run normalization
echo "2. Running normalization..."
echo "Start time: $(date)"

python "${PYTHON_SCRIPT}" \
    --ase-count-dir "${ASE_COUNT_DIR}" \
    --sample-list "${SAMPLE_LIST}" \
    --mapped-reads "${MAPPED_READS}" \
    --samplesheet "${SAMPLESHEET}" \
    --output "${OUT_TSV}" \
    --min-depth "${MIN_DEPTH}"

echo ""
echo "Normalization complete. Output saved to: ${OUT_TSV}"
echo "End time: $(date)"
