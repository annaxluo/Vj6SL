#!/bin/bash
#SBATCH --job-name=t5-1
#SBATCH --array=0-1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=t5-1_dexseq_count_%A_%a.out
#SBATCH --error=t5-1_dexseq_count_%A_%a.err

# ----------------------------------------
# Count reads aligned to DEXSeq exon bins. 
# 
# Usage:
#   sbatch t5-2_dexseq_count_mouse.sh
#   sbatch t5-2_dexseq_count_mouse.sh 2
# ----------------------------------------

export PYTHONUNBUFFERED=1

set -euo pipefail

# validate inputs
SAMPLES_PER_TASK="${1:-2}"
if ! [[ "${SAMPLES_PER_TASK}" =~ ^[0-9]+$ ]]; then
    echo "Error: SAMPLES_PER_TASK must be an integer"
    exit 1
fi

# load modules
module load conda/3-23.3.1
conda activate rnaseq
module load samtools/1.18

# paths 
BASE_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
PYTHON_SCRIPT_DIR="${BASE_DIR}/python/DEXSeq"
DEXSEQ_COUNT="${PYTHON_SCRIPT_DIR}/dexseq_count.py"

DEXSEQ_DIR="${BASE_DIR}/data/5_DEXSeq"
ANNOT_DIR="${DEXSEQ_DIR}/annotation"
COUNT_DIR="${DEXSEQ_DIR}/counts"
LOG_DIR="${DEXSEQ_DIR}/logs"

DEXSEQ_GFF="${ANNOT_DIR}/GRCm39.dexseq.gff"

# alignment paths 
ALIGN_DIR="${BASE_DIR}/data/3_STAR_Alignment"
FILTERED_BAM_DIR="${ALIGN_DIR}/filtered_bams"

SAMPLE_LIST="${BASE_DIR}/data/sample_list.txt"

# DEXSeq parameter settings
PAIRED="yes"
STRANDEDNESS="no"
SORT_ORDER="pos"

# output dirs
mkdir -p "${COUNT_DIR}" "${LOG_DIR}" 

ARRAY_TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"

# script info 
SCRIPT="t5-1_dexseq_count_mouse.sh"
echo "******************************"
echo "Script name: ${SCRIPT}"
echo "SLURM job ID: ${SLURM_JOB_ID:-local}"
echo "Start time: $(date)"
echo "******************************"
echo ""

# check files 
if [[ ! -f "${DEXSEQ_COUNT}" ]]; then
    echo "Error: Cannot find dexseq_count.py:"
    exit 1
fi

if [[ ! -f "${DEXSEQ_GFF}" ]]; then
    echo "Error: Cannot find DEXSeq annotation GFF:"
    exit 1
fi

# perform DEXSeq count for samples in each job array 
# task parameters
THREADS="${SLURM_CPUS_PER_TASK:-4}"
ARRAY_TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"

# check sample list 
if [[ ! -f "${SAMPLE_LIST}" ]]; then
    echo "Error: Sample list not found: ${SAMPLE_LIST}" 
    exit 1
fi

mapfile -t SAMPLES < "${SAMPLE_LIST}"
N_SAMPLES=$(wc -l < "${SAMPLE_LIST}")

START_INDEX=$(( ARRAY_TASK_ID * SAMPLES_PER_TASK ))
END_INDEX=$(( START_INDEX + SAMPLES_PER_TASK - 1 ))

if [[ "${START_INDEX}" -gt $(( N_SAMPLES - 1)) ]]; then
    echo "No samples assigned to this array task."
    echo "START_INDEX=${START_INDEX}, TOTAL_SAMPLES=${N_SAMPLES}"
    exit 0
fi

if [[ "${END_INDEX}" -gt $(( N_SAMPLES - 1)) ]]; then
    END_INDEX=$(( N_SAMPLES - 1))
fi

N_SAMPLES_TASK=$(( END_INDEX - START_INDEX + 1 ))

echo "Array task ID: ${ARRAY_TASK_ID}"
echo "Total samples: ${N_SAMPLES}"
echo "Number of samples in the task: ${N_SAMPLES_TASK}"
echo "Processing sample list lines: ${START_INDEX}-${END_INDEX}"
echo "-----------------------------"
echo ""

# Run DEXSeq counting for samples in each batch 
echo "Running DEXSeq counting..."
echo ""

process_sample(){

    local SAMPLE="$1"

    local BAM="${FILTERED_BAM_DIR}/${SAMPLE}.wasp_filtered.rg.bam"
    local OUT_COUNTS="${COUNT_DIR}/${SAMPLE}.dexseq_counts.txt"
    local OUT_LOG="${LOG_DIR}/${SAMPLE}.dexseq_count.log"
	
	echo "  Start time: $(date)"
	echo "  Input BAM file: ${BAM}"
	
	if [[ ! -f "${BAM}" ]]; then
        echo "BAM not found for sample ${SAMPLE}. Skipping."
        continue
    fi

    if [[ ! -f "${BAM}.bai" ]]; then
        echo "BAM index not found. Creating index:"
        echo "${BAM}.bai"
        samtools index "${BAM}"
    fi

    python "${DEXSEQ_COUNT}" \
        -p "${PAIRED}" \
        -s "${STRANDEDNESS}" \
        -r "${SORT_ORDER}" \
        "${DEXSEQ_GFF}" \
        "${BAM}" \
        "${OUT_COUNTS}" \
        > "${OUT_LOG}" 2>&1
	
	if [[ ! -s "${OUT_COUNTS}" ]]; then
        echo "Error: Unable to perform DEXSeq count for sample: ${SAMPLE}" 
        return 1
    fi
	
	echo "  Completed sample: ${SAMPLE}"
	echo "  Output saved to: ${OUT_COUNTS}"
	echo "  Log saved to: ${OUT_LOG}"
	echo "  End time: $(date)"
	echo "-----------------------------" 
	echo ""
	
	return 0
}

SAMPLE_NUM=0
FAILED_SAMPLES=()

for LINE_NUM in $(seq "${START_INDEX}" "${END_INDEX}"); do

    SAMPLE="${SAMPLES[$LINE_NUM]}"
	SAMPLE_NUM=$((SAMPLE_NUM + 1))

    echo "Processing sample ${SAMPLE_NUM}/${N_SAMPLES_TASK}: ${SAMPLE}"
	
    if ! process_sample "${SAMPLE}"; then
        echo "Error: sample failed: ${SAMPLE}" 
        FAILED_SAMPLES+=("${SAMPLE}")
    fi
done

# summary 
if [[ "${#FAILED_SAMPLES[@]}" -gt 0 ]]; then
    echo "The following samples failed:" 
    printf '%s\n' "${FAILED_SAMPLES[@]}" 
    exit 1
else
    echo "All samples in this array task completed successfully."
fi

echo ""
echo "Array task ${ARRAY_TASK_ID} complete"
echo "End time: $(date)"