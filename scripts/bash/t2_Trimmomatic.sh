#!/bin/bash
#SBATCH --job-name=t2
#SBATCH --output=t2_trim_%A_%a.out
#SBATCH --error=t2_trim_%A_%a.err
#SBATCH --array=0-1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00

# --------------------------------------
# Trim adapter from raw fastq files using Trimmomatic
# 1. Trim adapter
# 
# Usage:
#   sbatch t2_Trimmomatic.sh
#   sbatch t2_Trimmomatic.sh 2
# --------------------------------------

export MALLOC_ARENA_MAX=2

# validate inputs
SAMPLES_PER_TASK="${1:-2}"
if ! [[ "${SAMPLES_PER_TASK}" =~ ^[0-9]+$ ]]; then
    echo "Error: SAMPLES_PER_TASK must be an integer"
    exit 1
fi

# Define paths
BASE_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
INPUT_DIR="${BASE_DIR}/data/0_RawData"
OUTPUT_DIR="${BASE_DIR}/data/2_TrimmedData"
SAMPLE_LIST="${BASE_DIR}/data/sample_list.txt"

TRIMMOMATIC_JAR="${BASE_DIR}/external/Trimmomatic-0.39/trimmomatic-0.39.jar"
ADAPTERS="${BASE_DIR}/external/Trimmomatic-0.39/adapters/TruSeq3-PE.fa"

mkdir -p ${OUTPUT_DIR}

# script info 
SCRIPT="t2_Trimmomatic.sh"
echo "******************************"
echo "Script name: ${SCRIPT}"
echo "SLURM job ID: ${SLURM_JOB_ID:-local}"
echo "Start date and time: $(date)"
echo "******************************"
echo ""

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

# Process samples
echo "Running Trimmomatic..."
echo "Start time: $(date)"

SAMPLE_NUM=0
FAILED=0

for IDX in $(seq ${START_INDEX} ${END_INDEX}); do

    SAMPLE=${SAMPLES[$IDX]}
    SAMPLE_NUM=$((SAMPLE_NUM + 1))
    
    echo ""
    echo "Processing sample ${SAMPLE_NUM}/${N_SAMPLES_TASK}: ${SAMPLE}"

    R1_INPUT="${INPUT_DIR}/${SAMPLE}_1.fq.gz"
    R2_INPUT="${INPUT_DIR}/${SAMPLE}_2.fq.gz"

    # output files
    R1_PAIRED="${OUTPUT_DIR}/${SAMPLE}_1_paired.fq.gz"
    R1_UNPAIRED="${OUTPUT_DIR}/${SAMPLE}_1_unpaired.fq.gz"
    R2_PAIRED="${OUTPUT_DIR}/${SAMPLE}_2_paired.fq.gz"
    R2_UNPAIRED="${OUTPUT_DIR}/${SAMPLE}_2_unpaired.fq.gz"

    # Trimmomatic with default parameters for paired-end data
    java -XX:ActiveProcessorCount=${SLURM_CPUS_PER_TASK} \
        -jar ${TRIMMOMATIC_JAR} PE \
        -threads ${THREADS} \
        ${R1_INPUT} ${R2_INPUT} \
        ${R1_PAIRED} ${R1_UNPAIRED} \
        ${R2_PAIRED} ${R2_UNPAIRED} \
        ILLUMINACLIP:${ADAPTERS}:2:30:10 \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:4:15 \
        MINLEN:36

    # check exit status
    if [[ $? -eq 0 ]]; then
        echo "OK: ${SAMPLE}"
    else
        echo "Error: ${SAMPLE}"
        FAILED=$((FAILED + 1))
    fi

    echo "Sample complete time: $(date)"
	echo ""
done

# Print summary for this task
echo "----------------------------------"
echo "Task ${ARRAY_TASK_ID} complete"
echo "Samples processed: ${TASK_TOTAL}"
echo "Failed samples: ${FAILED}"
echo "End time: $(date)"
