#!/bin/bash
#SBATCH --job-name=t3-3
#SBATCH --array=0-1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=t3-3_add_rg_%A_%a.out
#SBATCH --error=t3-3_add_rg_%A_%a.err

# ----------------------------------------
# Add read groups to WASP-filtered BAMs  
# 1. Add read groups. 
# 
# Usage:
#   sbatch t3-3_ase_prepare.sh 
#   sbatch t3-3_ase_prepare.sh 2
# ----------------------------------------

set -euo pipefail

# validate inputs
SAMPLES_PER_TASK="${1:-2}"
if ! [[ "${SAMPLES_PER_TASK}" =~ ^[0-9]+$ ]]; then
    echo "Error: SAMPLES_PER_TASK must be an integer"
    exit 1
fi

set +u
module load samtools/1.18
module load gatk/4.5.0.0
set -u

# paths 
BASE_DIR="${SLURM_SUBMIT_DIR}"

ALIGN_DIR="${BASE_DIR}/data/3_STAR_Alignment"
FILTERED_BAM_DIR="${ALIGN_DIR}/filtered_bams"
STATS_DIR="${ALIGN_DIR}/stats"
REPORT_DIR="${ALIGN_DIR}/reports"

SAMPLE_LIST="${BASE_DIR}/data/sample_list.txt"

mkdir -p "${FILTERED_BAM_DIR}" "${STATS_DIR}" "${REPORT_DIR}"

# script info 
SCRIPT="t3-3_add_readGroups.sh"
echo "******************************"
echo "Script name: ${SCRIPT}"
echo "SLURM job ID: ${SLURM_JOB_ID:-local}"
echo "Start date and time: $(date)"
echo "******************************"
echo ""

# process each job array 
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

# add read group to one sample
process_sample(){

    local SAMPLE="$1"

    local IN_BAM="${FILTERED_BAM_DIR}/${SAMPLE}.wasp_filtered.bam"
    local OUT_BAM="${FILTERED_BAM_DIR}/${SAMPLE}.wasp_filtered.rg.bam"
    local REPORT="${REPORT_DIR}/${SAMPLE}.add_read_groups.report.txt"
    local FLAGSTAT="${STATS_DIR}/${SAMPLE}.wasp_filtered.rg.flagstat.txt"

    echo "  Input BAM: ${IN_BAM}"
    echo "  Output BAM: ${OUT_BAM}"

    if [[ ! -f "${IN_BAM}" ]]; then
        echo "  Error: input BAM not found for sample ${SAMPLE}: ${IN_BAM}" 
        return 1
    fi

    if [[ -f "${OUT_BAM}" && -f "${OUT_BAM}.bai" ]]; then
        echo "  Output BAM and index already exist. Skipping sample: ${SAMPLE}"
        echo "  Existing output: ${OUT_BAM}"
        return 0
    fi

    rm -f "${OUT_BAM}" "${OUT_BAM}.bai"

    gatk --java-options "-Xmx12g" AddOrReplaceReadGroups \
        -I "${IN_BAM}" \
        -O "${OUT_BAM}" \
        -RGID "${SAMPLE}" \
        -RGLB "bulk_RNAseq" \
        -RGPL "ILLUMINA" \
        -RGPU "${SAMPLE}" \
        -RGSM "${SAMPLE}" \
        --CREATE_INDEX false \
        --VALIDATION_STRINGENCY LENIENT

    samtools index -@ "${THREADS}" "${OUT_BAM}"

    # verify RG line
    samtools view -H "${OUT_BAM}" | grep '^@RG' || {
        echo "  Error: no @RG line found in output BAM for sample ${SAMPLE}" 
        return 1
    }

    samtools flagstat -@ "${THREADS}" "${OUT_BAM}" > "${FLAGSTAT}"

    {
        echo "Input BAM: ${IN_BAM}"
        echo "Output BAM: ${OUT_BAM}"
        echo
        echo "Read group added:"
        echo "  RGID=${SAMPLE}"
        echo "  RGLB=bulk_RNAseq"
        echo "  RGPL=ILLUMINA"
        echo "  RGPU=${SAMPLE}"
        echo "  RGSM=${SAMPLE}"
        echo
        echo "Output @RG header:"
        samtools view -H "${OUT_BAM}" | grep '^@RG' || true
        echo
        echo "Flagstat:"
        cat "${FLAGSTAT}"
    } > "${REPORT}"

    echo "  Finished sample: ${SAMPLE}"
    echo "  Output BAM: ${OUT_BAM}"
    echo "  Output index: ${OUT_BAM}.bai"
	echo "  End time: $(date)"
	echo "-------------------------------------"
	echo ""

    return 0
}

# process all samples in an array 
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
echo ""
echo "Array task ${ARRAY_TASK_ID} complete"
echo "End time: $(date)"

if [[ "${#FAILED_SAMPLES[@]}" -gt 0 ]]; then
    echo "The following samples failed:" 
    printf '%s\n' "${FAILED_SAMPLES[@]}" 
    exit 1
else
    echo "All samples in this array task completed successfully."
fi

