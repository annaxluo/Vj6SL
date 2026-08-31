#!/bin/bash
#SBATCH --job-name=t4-2
#SBATCH --output=t4-2_extract_variant_reads_%A_%a.out
#SBATCH --error=t4-2_extract_variant_reads_%A_%a.err
#SBATCH --array=0-1
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# ----------------------------------------
# Extract reads that align to the variant positions specified in a VCF file.  
# 1. Extract reads. 
# 
# Usage:
#   sbatch t4-2_ase_extract_variant_reads.sh
#   sbatch t4-2_ase_extract_variant_reads.sh 2
# ----------------------------------------

set -eo pipefail

# validate inputs
SAMPLES_PER_TASK="${1:-2}"

if ! [[ "${SAMPLES_PER_TASK}" =~ ^[0-9]+$ ]]; then
    echo "Error: SAMPLES_PER_TASK must be an integer"
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
OUT_DIR="${ASE_DIR}/variant_read_sequences"
VCF_FILE="${BASE_DIR}/variant_files/variant_list_mm39.vcf"

SAMPLE_LIST="${BASE_DIR}/data/sample_list.txt"

PY_SCRIPT="${BASE_DIR}/python/t4_extract_variant_reads.py"

mkdir -p "${OUT_DIR}"

# script info 
SCRIPT="t4-2_ase_extract_variant_reads.sh"
echo "******************************"
echo "Script name: ${SCRIPT}"
echo "SLURM job ID: ${SLURM_JOB_ID:-local}"
echo "Start date and time: $(date)"
echo "******************************"
echo ""

THREADS="${SLURM_CPUS_PER_TASK:-1}"
ARRAY_TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"

# check files 
if [[ ! -f "${PY_SCRIPT}" ]]; then
    echo "Error: Python script not found:"
    echo "  ${PY_SCRIPT}"
    exit 1
fi

if [[ ! -f "${VCF_FILE}" ]]; then
    echo "Error: VCF file not found:" 
    echo "  ${VCF_FILE}" 
    exit 1
fi

# read variants from VCF files 
mapfile -t VARIANT_LINES < <(
    awk 'BEGIN{OFS="\t"}
        /^#/ {next}
        {
            print $1, $2, $3, $4, $5
        }
    ' "${VCF_FILE}"
)

VARIANT_COUNT="${#VARIANT_LINES[@]}"

echo "number of variants from VCF: ${VARIANT_COUNT}"

# process each samples in each batch 
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

# count reads
SAMPLE_NUM=0
FAILED_SAMPLES=()

for LINE_NUM in $(seq "${START_INDEX}" "${END_INDEX}"); do

    SAMPLE="${SAMPLES[$LINE_NUM]}"
	SAMPLE_NUM=$((SAMPLE_NUM + 1))

    BAM="${FILTERED_BAM_DIR}/${SAMPLE}.wasp_filtered.rg.bam"
    BAI="${BAM}.bai"

    echo "Processing sample ${SAMPLE_NUM}/${N_SAMPLES_TASK}: ${SAMPLE}"
	echo "  Start time: $(date)"
    echo "  BAM: ${BAM}"

    # check input files
    if [[ ! -f "${BAM}" ]]; then
        echo "  Error: BAM not found for sample ${SAMPLE}. Skipping."
        FAILED_SAMPLES+=("${SAMPLE}")
        continue
    fi

    if [[ ! -f "${BAI}" ]]; then
        echo "  Error: BAM index not found for sample ${SAMPLE}. Skipping."
        FAILED_SAMPLES+=("${SAMPLE}")
        continue
    fi
	
    # extract reads for each variant
    for VARIANT_LINE in "${VARIANT_LINES[@]}"; do

        IFS=$'\t' read -r CHROM POS ID REF ALT <<< "${VARIANT_LINE}"

        SAFE_ID="${ID//[^A-Za-z0-9_.-]/_}"
        PREFIX="${OUT_DIR}/${SAMPLE}.${CHROM}_${POS}.${SAFE_ID}.${REF}_${ALT}"

        python "${PY_SCRIPT}" \
            "${BAM}" \
            "${CHROM}" \
            "${POS}" \
            "${REF}" \
            "${ALT}" \
            "${PREFIX}"

    done
	
	echo "  Complete."
	echo "-------------------------------"
	echo ""

done

# summary
if [[ "${#FAILED_SAMPLES[@]}" -gt 0 ]]; then
    echo ""
    echo "The following samples failed or were skipped:" 
    printf '%s\n' "${FAILED_SAMPLES[@]}" 
    exit 1
else
    echo "All samples in this array task completed successfully."
fi

echo "Array task ${ARRAY_TASK_ID} complete"
echo "End time: $(date)"