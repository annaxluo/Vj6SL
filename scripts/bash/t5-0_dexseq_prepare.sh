#!/bin/bash
#SBATCH --job-name=t5-0
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=t5-0_dexseq_prepare_%j.out
#SBATCH --error=t5-0_dexseq_prepare_%j.err

# ----------------------------------------
# Prepare gene annotations for DEXSeq exon read counts. 
# 
# Usage:
#   sbatch t5-1_dexseq_prepare.sh
# ----------------------------------------
export PYTHONUNBUFFERED=1

set -euo pipefail

module load conda/3-23.3.1
conda activate rnaseq

# paths 
BASE_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
PYTHON_SCRIPT_DIR="${BASE_DIR}/python"
DEXSEQ_PY_DIR="${PYTHON_SCRIPT_DIR}/DEXSeq" # update paths
PREPARE_ANNOT="${DEXSEQ_PY_DIR}/dexseq_prepare_annotation.py"

GTF="${BASE_DIR}/ref_data/refdata-gex-GRCm39-2024-A/genes/genes.gtf"

OUT_DIR="${BASE_DIR}/data/5_DEXSeq"
ANNOT_DIR="${OUT_DIR}/annotation"
DEXSEQ_GFF="${ANNOT_DIR}/GRCm39.dexseq.gff"

mkdir -p "${ANNOT_DIR}"

# script info 
SCRIPT="t5-0_dexseq_prepare.sh"
echo "******************************"
echo "Script name: ${SCRIPT}"
echo "SLURM job ID: ${SLURM_JOB_ID:-local}"
echo "Start date and time: $(date)"
echo "******************************"
echo ""

# check files 
if [[ ! -f "${PREPARE_ANNOT}" ]]; then
    echo "Error: Cannot find dexseq_prepare_annotation.py:"
    exit 1
fi

if [[ ! -f "${GTF}" ]]; then
    echo "Error: Cannot find GTF:"
    exit 1
fi

# Run DEXSeq annotation preparation
echo "Preparing DEXSeq annotation..."
echo "Start time: $(date)"
echo "  Input GTF:  ${GTF}"

python "${PREPARE_ANNOT}" \
    "${GTF}" \
    "${DEXSEQ_GFF}"

echo "  DEXSeq annotation complete. Annotations saved to: ${DEXSEQ_GFF}"
echo "End time: $(date)"
