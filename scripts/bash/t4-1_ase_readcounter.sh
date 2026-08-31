#!/bin/bash
#SBATCH --job-name=t4-1
#SBATCH --array=0-1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=04:00:00
#SBATCH --output=t4-1_ase_readcounter_%A_%a.out
#SBATCH --error=t4-1_ase_readcounter_%A_%a.err

# ----------------------------------------
# Count reads for each allele.   
# 1. Verify required input files.  
# 2. Count allele expression using ASEReadCounter. 
# 
# Usage:
#   sbatch t4-1_ase_readcounter.sh 
#   sbatch t4-1_ase_readcounter.sh 2
#   sbatch t4-1_ase_readCounter.sh 2 chr12
#   sbatch t4-1_ase_readCounter.sh 2 chr12:118260000-118270000
# ----------------------------------------

set -euo pipefail

# validate inputs
SAMPLES_PER_TASK="${1:-2}"
ASE_INTERVAL="${2:-chr12}"

if ! [[ "${SAMPLES_PER_TASK}" =~ ^[0-9]+$ ]]; then
    echo "Error: SAMPLES_PER_TASK must be an integer"
    exit 1
fi

set +u
module load samtools/1.18
module load gatk/4.5.0.0
module load htslib/1.18 
module load bcftools/1.18
set -u

# paths 
BASE_DIR="${SLURM_SUBMIT_DIR:-$PWD}"

ALIGN_DIR="${BASE_DIR}/data/3_STAR_Alignment"
FILTERED_BAM_DIR="${ALIGN_DIR}/filtered_bams"

ASE_DIR="${BASE_DIR}/data/4_ASE_Analysis"
COUNT_DIR="${ASE_DIR}/allele_counts_Q20"
REPORT_DIR="${ASE_DIR}/reports"
STATS_DIR="${ASE_DIR}/stats"

SAMPLE_LIST="${BASE_DIR}/data/sample_list.txt"

# reference data, variant data 
REF_FASTA="${BASE_DIR}/ref_data/refdata-gex-GRCm39-2024-A/fasta/genome.fa"
REF_DICT="${REF_FASTA%.*}.dict"
VCF_FILE="${BASE_DIR}/variant_files/variant_list_mm39.sorted.vcf.gz"

mkdir -p "${COUNT_DIR}" "${REPORT_DIR}" "${STATS_DIR}"

# script info 
SCRIPT="t4-1_ase_readcounter.sh"
echo "******************************"
echo "Script name: ${SCRIPT}"
echo "SLURM job ID: ${SLURM_JOB_ID:-local}"
echo "Start date and time: $(date)"
echo "******************************"
echo ""

# 1. verify inputs
if [[ ! -f "${SAMPLE_LIST}" ]]; then
    echo "Error: sample list not found: ${SAMPLE_LIST}" 
    exit 1
fi

if [[ ! -f "${VCF_FILE}" ]]; then
    echo "Error: VCF file not found: ${VCF_FILE}" 
    exit 1
fi

if [[ ! -f "${VCF_FILE}.tbi" ]]; then
    echo "Error: VCF index not found: ${VCF_FILE}.tbi" 
    exit 1
fi

if [[ ! -f "${REF_FASTA}" ]]; then
    echo "Error: reference FASTA not found: ${REF_FASTA}" 
    exit 1
fi

if [[ ! -f "${REF_FASTA}.fai" ]]; then
    echo "Reference FASTA index not found. ${REF_FASTA}.fai"
    exit 1
fi

if [[ ! -f "${REF_DICT}" ]]; then
    echo "Reference dictionary not found. ${REF_DICT}"
    exit 1
fi

# 1. Run ASEReadCounter for samples in each batch  
# task parameters
THREADS="${SLURM_CPUS_PER_TASK:-2}"
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

SAMPLE_NUM=0
FAILED_SAMPLES=()

for LINE_NUM in $(seq "${START_INDEX}" "${END_INDEX}"); do

    SAMPLE="${SAMPLES[$LINE_NUM]}"
	SAMPLE_NUM=$((SAMPLE_NUM + 1))

    BAM="${FILTERED_BAM_DIR}/${SAMPLE}.wasp_filtered.rg.bam" # with read groups
    OUT="${COUNT_DIR}/${SAMPLE}.ASEReadCounter.tsv"

    REPORT="${REPORT_DIR}/${SAMPLE}.ASEReadCounter.report.txt"

    echo "Processing sample ${SAMPLE_NUM}/${N_SAMPLES_TASK}: ${SAMPLE}"
	echo "  Start time: $(date)"
    echo "  BAM: ${BAM}"
    echo "  Output: ${OUT}"

    if [[ ! -f "${BAM}" ]]; then
        echo "  Error: BAM not found: ${BAM}" 
        FAILED_SAMPLES+=("${SAMPLE}")
        continue
    fi

    if [[ ! -f "${BAM}.bai" ]]; then
        echo "  BAM index not found. Creating index."
        samtools index -@ "${THREADS}" "${BAM}"
    fi

    # check read group
    samtools view -H "${BAM}" | grep -q '^@RG' || {
        echo "  Error: no @RG line found in BAM for ${SAMPLE}" 
        FAILED_SAMPLES+=("${SAMPLE}")
        continue
    }

    echo "  Running ASEReadCounter on interval: ${ASE_INTERVAL}"

    gatk --java-options "-Xmx10g" ASEReadCounter \
        -R "${REF_FASTA}" \
        -I "${BAM}" \
        -V "${VCF_FILE}" \
        -O "${OUT}" \
        -L "${ASE_INTERVAL}" \
        --min-depth 1 \
        --min-base-quality 20 \
        --min-mapping-quality 20

    echo "  ASEReadCounter completed."

	{
		echo "Sample: ${SAMPLE}"
		echo
		echo "BAM: ${BAM}"
		echo "VCF: ${VCF_FILE}"
		echo "Reference: ${REF_FASTA}"
		echo "ASEReadCounter interval: ${ASE_INTERVAL}"
		echo "Output: ${OUT}"
		echo
		echo "ASEReadCounter command:"
		echo "gatk ASEReadCounter \\"
		echo "  -R ${REF_FASTA} \\"
		echo "  -I ${BAM} \\"
		echo "  -V ${VCF_FILE} \\"
		echo "  -O ${OUT} \\"
		echo "  -L ${ASE_INTERVAL} \\"
		echo "  --min-depth 1 \\"
		echo "  --min-base-quality 20 \\"
		echo "  --min-mapping-quality 20"
    } > "${REPORT}"

    echo "  Done sample: ${SAMPLE}"
    echo "  Output: ${OUT}"
    echo "  Report: ${REPORT}"
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

echo "Array task ${ARRAY_TASK_ID} complete"
echo "End time: $(date)"
