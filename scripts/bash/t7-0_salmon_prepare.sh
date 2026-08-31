#!/bin/bash
#SBATCH --job-name=t7-0
#SBATCH --output=t7-0_salmon_prepare_%j.out
#SBATCH --error=t7-0_salmon_prepare_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G

# ----------------------------------------
# Build Salmon transcriptome index. 
# 1. Check input files 
# 2. Extract transcript sequences from GTF 
# 3. Build Salmon index
# 
# Usage:
#   sbatch t7-1_salmon_prepare.sh
# ----------------------------------------

set -euo pipefail

module load Salmon/1.10.1
module load gffread/0.12.7

# paths
BASE_DIR="${SLURM_SUBMIT_DIR:-$PWD}"

REF_DATA_DIR="${BASE_DIR}/ref_data/refdata-gex-GRCm39-2024-A"
GENOME_FA="${REF_DATA_DIR}/fasta/genome.fa"
GTF="${REF_DATA_DIR}/genes/genes.gtf"

OUT_DIR="${BASE_DIR}/data/7_Salmon"
REF_DIR="${OUT_DIR}/reference"
INDEX_DIR="${OUT_DIR}/index/GRCm39_salmon_index"
LOG_DIR="${OUT_DIR}/logs"

TRANSCRIPTS_FA="${REF_DIR}/GRCm39.transcripts.fa"

mkdir -p "${REF_DIR}" "${INDEX_DIR}" "${LOG_DIR}" 

# script info 
THREADS="${SLURM_CPUS_PER_TASK:-1}"
SCRIPT="t7-0_salmon_prepare.sh"

echo "******************************"
echo "Script name: ${SCRIPT}"
echo "SLURM job ID: ${SLURM_JOB_ID:-local}"
echo "Start time: $(date)"
echo "******************************"
echo ""

# 1. check input files
if [[ ! -f "${GENOME_FA}" ]]; then
    echo "Error: genome FASTA not found: ${GENOME_FA}"
    exit 1
fi

if [[ ! -f "${GTF}" ]]; then
    echo "Error: GTF not found: ${GTF}"
    exit 1
fi

# 2. Extract transcript FASTA
if [[ ! -f "${TRANSCRIPTS_FA}" ]]; then
    echo ""
    echo "Extracting transcript sequences..."

    gffread "${GTF}" \
        -g "${GENOME_FA}" \
        -w "${TRANSCRIPTS_FA}" \
        > "${LOG_DIR}/gffread_transcripts.log" \
        2> "${LOG_DIR}/gffread_transcripts.err"
else
    echo ""
    echo "Transcript FASTA exists:  ${TRANSCRIPTS_FA}"
fi
echo ""

# 3. Build Salmon index
echo "Building Salmon index..."

salmon index \
    -t "${TRANSCRIPTS_FA}" \
    -i "${INDEX_DIR}" \
    -p "${THREADS}" \
    --keepDuplicates \
    > "${LOG_DIR}/salmon_index.log" \
    2> "${LOG_DIR}/salmon_index.err"

echo "Salmon index saved to: ${INDEX_DIR}"
echo "End time: $(date)"
