#!/bin/bash
#SBATCH --job-name=t3-4
#SBATCH --array=0-1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=t3-4_picard_qc_%A_%a.out
#SBATCH --error=t3-4_picard_qc_%A_%a.err

# ----------------------------------------
# Perform QC using gatk/Picard tools on the WASP-filtered BAM files. 
# 1. ValidateSamFile
# 2. CollectAlignmentSummaryMetrics
# 3. CollectInsertSizeMetrics
# 4. MarkDuplicates metrics only
# 5. CollectRnaSeqMetrics
#
# Usage:
#   sbatch t3-4_picard_qc.sh
#   sbatch t3-4_picard_qc.sh 2
# ----------------------------------------

set -euo pipefail

# validate inputs
SAMPLES_PER_TASK="${1:-2}"
if ! [[ "${SAMPLES_PER_TASK}" =~ ^[0-9]+$ ]]; then
    echo "Error: SAMPLES_PER_TASK must be an integer"
    exit 1
fi

set +u
module load gatk/4.5.0.0
module load samtools/1.18
set -u

# paths
BASE_DIR="${SLURM_SUBMIT_DIR:-$PWD}"

ALIGN_DIR="${BASE_DIR}/data/3_STAR_Alignment"
FILTERED_BAM_DIR="${ALIGN_DIR}/filtered_bams"

SAMPLE_LIST="${BASE_DIR}/data/sample_list.txt"

OUT_DIR="${ALIGN_DIR}/multiqc_outputs"
PICARD_DIR="${OUT_DIR}/picard"
TMP_DIR="${OUT_DIR}/tmp"

STATS_DIR="${ALIGN_DIR}/stats"
REPORT_DIR="${ALIGN_DIR}/reports"

# Reference files
REF_FASTA="${BASE_DIR}/ref_data/refdata-gex-GRCm39-2024-A/fasta/genome.fa"
REF_FLAT="${BASE_DIR}/ref_data/refFlat/GRCm39/refFlat.txt"
# RNA-seq strandedness for Picard CollectRnaSeqMetrics
STRAND_SPECIFICITY="NONE"

mkdir -p "${PICARD_DIR}" "${TMP_DIR}" "${STATS_DIR}" "${REPORT_DIR}"

# script info 
SCRIPT="t3-4_picard_qc.sh"
echo "******************************"
echo "Script name: ${SCRIPT}"
echo "SLURM job ID: ${SLURM_JOB_ID}"
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

# function for Picard QC
run_picard_qc(){

    local BAM="$1"
    local SAMPLE="$2"
    local STAGE="$3"
    local STAGE_OUT="$4"
	
	local SAMPLE_OUT="${STAGE_OUT}/${SAMPLE}"
    local SAMPLE_TMP="${TMP_DIR}/${SAMPLE}.${STAGE}"
	
	mkdir -p "${SAMPLE_OUT}" "${SAMPLE_TMP}"

	echo "  Start time: $(date)"
    echo "  Stage: ${STAGE}"
    echo "  BAM: ${BAM}"
	
	mkdir -p "${STAGE_OUT}/${SAMPLE}"

    # check files 
    if [[ ! -f "${BAM}" ]]; then
        echo "BAM not found. Skipping:"
        echo "${BAM}"
        return 0
    fi

    if [[ ! -f "${BAM}.bai" ]]; then
        echo "BAM index not found. Creating index:"
        echo "${BAM}.bai"
        samtools index "${BAM}"
    fi 

    # 1. ValidateSamFile
    echo "  1. Running ValidateSamFile..."

    set +e

    gatk ValidateSamFile \
        -I "${BAM}" \
        -O "${SAMPLE_OUT}/${SAMPLE}.${STAGE}.ValidateSamFile.txt" \
        --MODE SUMMARY \
        --TMP_DIR "${SAMPLE_TMP}" \
        --VALIDATION_STRINGENCY LENIENT

    VALIDATE_EXIT_CODE=$?

    set -e

    if [[ ${VALIDATE_EXIT_CODE} -ne 0 ]]; then
        echo "    ValidateSamFile reported validation errors for ${SAMPLE}."
        echo "    Exit code: ${VALIDATE_EXIT_CODE}"
    fi

    # 2. CollectAlignmentSummaryMetrics
    echo "  2. Running CollectAlignmentSummaryMetrics..."

    gatk CollectAlignmentSummaryMetrics \
        -R "${REF_FASTA}" \
        -I "${BAM}" \
        -O "${SAMPLE_OUT}/${SAMPLE}.${STAGE}.alignment_summary_metrics.txt" \
        --TMP_DIR "${SAMPLE_TMP}" \
        --VALIDATION_STRINGENCY LENIENT

    # 3. CollectInsertSizeMetrics
    echo "  3. Running CollectInsertSizeMetrics..."

    gatk CollectInsertSizeMetrics \
        -I "${BAM}" \
        -O "${SAMPLE_OUT}/${SAMPLE}.${STAGE}.insert_size_metrics.txt" \
        -H "${SAMPLE_OUT}/${SAMPLE}.${STAGE}.insert_size_histogram.pdf" \
        -M 0.5 \
        --TMP_DIR "${SAMPLE_TMP}" \
        --VALIDATION_STRINGENCY LENIENT

    # 4. MarkDuplicates metrics only
    echo "  4. Running MarkDuplicates..."

    gatk MarkDuplicates \
        -I "${BAM}" \
        -O "${SAMPLE_OUT}/${SAMPLE}.${STAGE}.markdup.tmp.bam" \
        -M "${SAMPLE_OUT}/${SAMPLE}.${STAGE}.duplication_metrics.txt" \
        --CREATE_INDEX false \
        --REMOVE_DUPLICATES false \
        --ASSUME_SORTED true \
        --TMP_DIR "${SAMPLE_TMP}" \
        --VALIDATION_STRINGENCY LENIENT

    rm -f "${SAMPLE_OUT}/${SAMPLE}.${STAGE}.markdup.tmp.bam"
    rm -f "${SAMPLE_OUT}/${SAMPLE}.${STAGE}.markdup.tmp.bai"
    rm -f "${SAMPLE_OUT}/${SAMPLE}.${STAGE}.markdup.tmp.bam.bai"

    # 5. CollectRnaSeqMetrics
    echo "  5. Running CollectRnaSeqMetrics..."

    gatk CollectRnaSeqMetrics \
	    -I "${BAM}" \
		-O "${SAMPLE_OUT}/${SAMPLE}.${STAGE}.rna_seq_metrics.txt" \
		--REF_FLAT "${REF_FLAT}" \
		--STRAND_SPECIFICITY "${STRAND_SPECIFICITY}" \
		--TMP_DIR "${SAMPLE_TMP}" \
		--VALIDATION_STRINGENCY LENIENT

    rm -rf "${SAMPLE_TMP}" 
    echo "  Completed Picard QC for ${SAMPLE}"
	echo "  End time: $(date)"
	echo ""
	
	return 0
}

# process samples in each array 
SAMPLE_NUM=0
FAILED_SAMPLES=()

for LINE_NUM in $(seq "${START_INDEX}" "${END_INDEX}"); do
    
	SAMPLE="${SAMPLES[$LINE_NUM]}"
    SAMPLE_NUM=$(( SAMPLE_NUM+1 ))
	
    echo "Processing sample ${SAMPLE_NUM}/${N_SAMPLES_TASK}: ${SAMPLE}"
		
    BAM="${FILTERED_BAM_DIR}/${SAMPLE}.wasp_filtered.rg.bam"

    if ! run_picard_qc \
        "${BAM}" \
        "${SAMPLE}" \
        "post_wasp_rg" \
        "${PICARD_DIR}"
    then
        echo "Error: Picard QC failed for sample: ${SAMPLE}" 
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

echo "Picard QC completed for array task ${ARRAY_TASK_ID}"
