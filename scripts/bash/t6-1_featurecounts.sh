#!/bin/bash
#SBATCH --job-name=t6-1
#SBATCH --array=0-1
#SBATCH --output=t6-1_featurecounts_%A_%a.out
#SBATCH --error=t6-1_featurecounts_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G

# ----------------------------------------
# Use featureCounts to extract exon and gene read counts from WASP-filtered BAMs
# 1. Build BAM lists.
# 2. Run featureCounts to estimate gene-level counts.
# 3. Run featureCounts to estimate exon-level counts. 
# 
# Usage:
#   sbatch t6-1_featureCounts.sh
#   sbatch t6-1_featureCounts.sh 2
# ----------------------------------------

set -euo pipefail

# validate inputs
SAMPLES_PER_TASK="${1:-2}"
if ! [[ "${SAMPLES_PER_TASK}" =~ ^[0-9]+$ ]]; then
    echo "Error: SAMPLES_PER_TASK must be an integer"
    exit 1
fi

module load subread/2.0.0

# paths 
BASE_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
ALIGN_DIR="${BASE_DIR}/data/3_STAR_Alignment"
BAM_DIR="${ALIGN_DIR}/filtered_bams"

GTF="${BASE_DIR}/ref_data/refdata-gex-GRCm39-2024-A/genes/genes.gtf"

OUT_DIR="${BASE_DIR}/data/6_featureCounts"
GENE_DIR="${OUT_DIR}/gene_counts"
EXON_DIR="${OUT_DIR}/exon_counts"

LOG_DIR="${OUT_DIR}/logs"
BAM_LIST="${OUT_DIR}/all_samples_bams.txt"

SAMPLE_LIST="${BASE_DIR}/data/sample_list.txt"

mkdir -p "${GENE_DIR}" "${EXON_DIR}" "${LOG_DIR}"

# featureCounts parameters
STRANDEDNESS=0
PAIRED_OPTS="-p -C"
BAM_SUFFIX=".wasp_filtered.rg.bam"

# script info 
SCRIPT="t6-1_featureCounts.sh"
echo "******************************"
echo "Script name: ${SCRIPT}"
echo "SLURM job ID: ${SLURM_JOB_ID:-local}"
echo "Start time: $(date)"
echo "******************************"
echo ""

# perform FeatureCounts for samples in each job array 
# task parameters
THREADS="${SLURM_CPUS_PER_TASK:-8}"
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

# Run featureCounts for each sample
SAMPLE_NUM=0
FAILED_SAMPLES=()

for LINE_NUM in $(seq "${START_INDEX}" "${END_INDEX}"); do

    SAMPLE="${SAMPLES[$LINE_NUM]}"
	SAMPLE_NUM=$((SAMPLE_NUM + 1))

    echo "Processing sample ${SAMPLE_NUM}/${N_SAMPLES_TASK}: ${SAMPLE}"
	
    BAM="${BAM_DIR}/${SAMPLE}${BAM_SUFFIX}"

    GENE_OUT="${GENE_DIR}/${SAMPLE}.gene_counts.txt"
    EXON_OUT="${EXON_DIR}/${SAMPLE}.exon_counts.txt"

    GENE_LOG="${LOG_DIR}/${SAMPLE}.featureCounts_gene.log"
    GENE_ERR="${LOG_DIR}/${SAMPLE}.featureCounts_gene.err"
    EXON_LOG="${LOG_DIR}/${SAMPLE}.featureCounts_exon.log"
    EXON_ERR="${LOG_DIR}/${SAMPLE}.featureCounts_exon.err"
	
	if [[ ! -f "${BAM}" ]]; then
        echo "Error: BAM not found for sample ${SAMPLE}."
        FAILED_SAMPLES+=("${SAMPLE}")
        continue
    fi

    # 1. Run featureCounts: gene-level counts
	echo "  Running gene-level counts..."
	echo "    Start time: $(date)"
	
	if featureCounts \
        -T "${THREADS}" \
        ${PAIRED_OPTS} \
        -s "${STRANDEDNESS}" \
        -a "${GTF}" \
        -t exon \
        -g gene_id \
        -o "${GENE_OUT}" \
        "${BAM}" \
        > "${GENE_LOG}" \
        2> "${GENE_ERR}"
    then
        echo "    Gene-level featureCounts completed."
        echo "    Output saved to: ${GENE_OUT}"
        echo "    Summary saved to: ${GENE_OUT}.summary"
        
    else
        echo "    Error: gene-level featureCounts failed."
        FAILED_SAMPLES+=("${SAMPLE}")
        continue
    fi

    echo "    End time: $(date)"
	echo ""
	
	# 2. Run featureCounts: exon-level counts
	echo "  Running exon-level counts..."
    echo "    Start time: $(date)"
    
	if featureCounts \
        -T "${THREADS}" \
        ${PAIRED_OPTS} \
        -s "${STRANDEDNESS}" \
        -f \
        -a "${GTF}" \
        -t exon \
        -g gene_id \
        -o "${EXON_OUT}" \
        "${BAM}" \
        > "${EXON_LOG}" \
        2> "${EXON_ERR}"
    then
        echo "    Exon-level featureCounts completed."
        echo "    Output saved to: ${EXON_OUT}"
        echo "    Summary saved to: ${EXON_OUT}.summary"
    else
        echo "    Error: exon-level featureCounts failed."
        FAILED_SAMPLES+=("${SAMPLE}")
        continue
    fi

    echo "    End time: $(date)"
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