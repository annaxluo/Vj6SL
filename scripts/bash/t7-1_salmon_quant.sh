#!/bin/bash
#SBATCH --job-name=t7-1
#SBATCH --array=0-1
#SBATCH --output=t7-1_salmon_quant_%A_%a.out
#SBATCH --error=t7-1_salmon_quant_%A_%a.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G

# ----------------------------------------
# Quantify transcript abundance with Salmon. 
# 
# 1. Check input files 
# 2. Extract transcript sequences from GTF 
# 3. Build Salmon index
# 
# Usage:
#   sbatch t7-2_salmon_quant.sh
#   sbatch t7-2_salmon_quant.sh 2
# ----------------------------------------

set -euo pipefail

# validate inputs
SAMPLES_PER_TASK="${1:-2}"
if ! [[ "${SAMPLES_PER_TASK}" =~ ^[0-9]+$ ]]; then
    echo "Error: SAMPLES_PER_TASK must be an integer"
    exit 1
fi

module load Salmon/1.10.1

# paths
BASE_DIR="${SLURM_SUBMIT_DIR:-$PWD}"

FASTQ_DIR="${BASE_DIR}/data/2_TrimmedData"

SAMPLE_LIST="${BASE_DIR}/data/sample_list.txt"

SALMON_DIR="${BASE_DIR}/data/7_Salmon"
INDEX_DIR="${SALMON_DIR}/index/GRCm39_salmon_index"
QUANT_DIR="${SALMON_DIR}/quant"
LOG_DIR="${SALMON_DIR}/logs"

mkdir -p "${QUANT_DIR}" "${LOG_DIR}"

# script info 
SCRIPT="t7-1_salmon_quant.sh"
echo "******************************"
echo "Script name: ${SCRIPT}"
echo "SLURM job ID: ${SLURM_JOB_ID:-local}"
echo "Start time: $(date)"
echo "******************************"
echo ""

ARRAY_TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
THREADS="${SLURM_CPUS_PER_TASK:-1}"

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

# Salmon quantification for each sample
SAMPLE_NUM=0
FAILED_SAMPLES=()

for LINE_NUM in $(seq "${START_INDEX}" "${END_INDEX}"); do

    SAMPLE="${SAMPLES[$LINE_NUM]}"
	SAMPLE_NUM=$((SAMPLE_NUM + 1))
	
	R1="${FASTQ_DIR}/${SAMPLE}_1_paired.fq.gz"
    R2="${FASTQ_DIR}/${SAMPLE}_2_paired.fq.gz"
    SAMPLE_QUANT_DIR="${QUANT_DIR}/${SAMPLE}"
	
	echo "Processing sample ${SAMPLE_NUM}/${N_SAMPLES_TASK}: ${SAMPLE}"
	echo "  Start time: $(date)"
	echo "  R1: ${R1}"
	echo "  R2: ${R2}"
	
	if [[ ! -f "${R1}" ]]; then
        echo "Error: R1 fastq not found: ${R1}"
        FAILED_SAMPLES+=("${SAMPLE}")
		continue
    fi

    if [[ ! -f "${R2}" ]]; then
        echo "Error: R2 fastq not found: ${R2}"
        FAILED_SAMPLES+=("${SAMPLE}")
		continue
    fi

    mkdir -p "${SAMPLE_QUANT_DIR}"
	
	set +e
    salmon quant \
        -i "${INDEX_DIR}" \
        -l A \
        -1 "${R1}" \
        -2 "${R2}" \
        -p "${THREADS}" \
        --validateMappings \
        --gcBias \
        --seqBias \
        -o "${SAMPLE_QUANT_DIR}" \
        > "${LOG_DIR}/${SAMPLE}.salmon_quant.log" \
        2> "${LOG_DIR}/${SAMPLE}.salmon_quant.err"
    
	EXIT_CODE=$? 
	
	set -e
	
	if [[ "${EXIT_CODE}" -ne 0 ]]; then
        echo "  Error: Salmon quant failed for sample ${SAMPLE}." 
        echo "  Exit code: ${EXIT_CODE}" 
        FAILED_SAMPLES+=("${SAMPLE}")
        continue
    fi

    echo "  Completed sample: ${SAMPLE}"
	echo "  End time: $(date)"
    echo "-----------------------------"
    echo "" 	
	
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
