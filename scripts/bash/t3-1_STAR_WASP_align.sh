#!/bin/bash
#SBATCH --job-name=t3-1
#SBATCH --output=t3-1_star_wasp_%A_%a.out
#SBATCH --error=t3-1_star_wasp_%A_%a.err
#SBATCH --array=0-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=24:00:00

# ----------------------------------------
# Perform alignment using STAR in WASP mode  
# 1. aligned trimmed fastq to reference
# 2. index output BAM files 
# 
# Usage:
#   sbatch t3-1_STAR_align.sh
#   sbatch t3-1_STAR_align.sh 2
# ----------------------------------------

set -euo pipefail

# validate inputs
SAMPLES_PER_TASK="${1:-2}"
if ! [[ "${SAMPLES_PER_TASK}" =~ ^[0-9]+$ ]]; then
    echo "Error: SAMPLES_PER_TASK must be an integer"
    exit 1
fi

# load modules
module load star/2.7.8a
module load samtools/1.18

# paths and variables
BASE_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
FASTQ_DIR="${BASE_DIR}/data/2_TrimmedData"
STAR_INDEX="${BASE_DIR}/ref_data/refdata-gex-GRCm39-2024-A/star_2.7.8a"
VCF_FILE="${BASE_DIR}/variant_files/variant_list_mm39.vcf"
OUTPUT_DIR="${BASE_DIR}/data/3_STAR_Alignment/aligned_bams"

SAMPLE_LIST="${BASE_DIR}/data/sample_list.txt"

mkdir -p "${STAR_INDEX_DIR}" "${OUTPUT_DIR}"

# script info 
SCRIPT="t3-1_STAR_WASP_align.sh"
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

process_sample(){

    local SAMPLE="$1"

    local R1="${FASTQ_DIR}/${SAMPLE}_1_paired.fq.gz"
    local R2="${FASTQ_DIR}/${SAMPLE}_2_paired.fq.gz"
    local SAMPLE_OUT_DIR="${OUTPUT_DIR}/${SAMPLE}"
	local BAM="${SAMPLE_OUT_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam"
	
	mkdir -p "${SAMPLE_OUT_DIR}"
	
	echo "  Running STAR alignment..."
	
	STAR --runThreadN ${THREADS} \
         --genomeDir ${STAR_INDEX} \
         --readFilesIn ${R1} ${R2} \
         --readFilesCommand zcat \
         --outFileNamePrefix ${SAMPLE_OUT_DIR}/${SAMPLE}_ \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMunmapped Within \
         --outSAMattributes NH HI AS nM NM MD jM jI rB MC vA vG vW \
         --alignEndsType EndToEnd \
         --waspOutputMode SAMtag \
         --varVCFfile ${VCF_FILE} \
         --outFilterMultimapNmax 1 \
         --alignIntronMin 20 \
         --alignIntronMax 1000000 \
         --alignMatesGapMax 1000000 \
         --limitBAMsortRAM 40000000000
		 
    # check if STAR completed successfully
	if [[ $? -eq 0 && -f "${SAMPLE_OUT_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam" ]]; then
        echo "  Alignment completed successfully"
        
		# index BAM file
        echo "  Indexing: ${BAM}"
        samtools index -@ 4 "${BAM}"
        
        if [[ -f "${BAM}.bai" ]]; then
            echo "  Index created: ${BAM}.bai"
        else
            echo "  Failed to create index for ${SAMPLE}"
        fi
		
    else
        echo "  Error: Alignment failed for ${SAMPLE}"
		return 1
    fi
	
	echo "  End time: $(date)"
    echo "---------------------------------"
	echo ""
	
	return 0

}

# process samples in the array 
SAMPLE_NUM=0
FAILED_SAMPLES=()

for LINE_NUM in $(seq "${START_INDEX}" "${END_INDEX}"); do
    
	SAMPLE="${SAMPLES[$LINE_NUM]}"
	SAMPLE_NUM=$((SAMPLE_NUM + 1))
	
    echo "Processing sample ${SAMPLE_NUM}/${N_SAMPLES_TASK}: ${SAMPLE}"
    echo "  Start time: $(date)"

    if ! process_sample "${SAMPLE}"; then
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