#!/bin/bash
#SBATCH --job-name=t3-2
#SBATCH --array=0-1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=t3-2_wasp_filter_%A_%a.out
#SBATCH --error=t3-2_wasp_filter_%A_%a.err

# ----------------------------------------
# Filter alignments by WASP tags 
# 1. filter alignments to retain alignments with WASP vW=1 or no vW tags
# 2. index output BAM files 
# 
# Usage:
#   sbatch t3-2_WASP_filter.sh
#   sbatch t3-2_WASP_filter.sh 2
# ----------------------------------------

set -euo pipefail

# validate inputs
SAMPLES_PER_TASK="${1:-2}"
if ! [[ "${SAMPLES_PER_TASK}" =~ ^[0-9]+$ ]]; then
    echo "Error: SAMPLES_PER_TASK must be an integer"
    exit 1
fi

# load modules
module load samtools/1.18

BASE_DIR="${SLURM_SUBMIT_DIR:-$PWD}"

ALIGN_DIR="${BASE_DIR}/data/3_STAR_Alignment"
INPUT_BAM_DIR="${ALIGN_DIR}/aligned_bams"
FILTERED_BAM_DIR="${ALIGN_DIR}/filtered_bams"
STATS_DIR="${ALIGN_DIR}/stats"
REPORT_DIR="${ALIGN_DIR}/reports"

SAMPLE_LIST="${BASE_DIR}/data/sample_list.txt"

mkdir -p "${FILTERED_BAM_DIR}" "${STATS_DIR}" "${REPORT_DIR}"

# script info 
SCRIPT="t3-2_WASP_filter.sh"
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

    local IN_BAM="${ALIGN_DIR}/${SAMPLE}/${SAMPLE}_Aligned.sortedByCoord.out.bam"
    local OUT_BAM="${FILTERED_BAM_DIR}/${SAMPLE}.wasp_filtered.bam"
    local TMP_BAM="${FILTERED_BAM_DIR}/${SAMPLE}.wasp_filtered.tmp.bam"
    local REPORT="${REPORT_DIR}/${SAMPLE}.wasp_filter_report.txt"

    echo "  Start time: $(date)"
    echo "  Input BAM: ${IN_BAM}"
    echo "  Output BAM: ${OUT_BAM}"
	echo ""

    # check files
    if [[ ! -f "${IN_BAM}" ]]; then
        echo "Error: Input BAM not found: ${IN_BAM}" 
        return 1
    fi

    if [[ ! -f "${IN_BAM}.bai" ]]; then
        echo "Input BAM index not found. Creating index for input BAM."
        samtools index -@ "${THREADS}" "${IN_BAM}"
    fi

    # count reads
    echo "  Counting input alignments..."

    local TOTAL_IN
    local NO_VW_IN
    local VW1_IN
    local VW_OTHER_IN

    TOTAL_IN=$(samtools view -@ "${THREADS}" -c "${IN_BAM}")

    NO_VW_IN=$(samtools view -@ "${THREADS}" "${IN_BAM}" | \
        awk 'BEGIN{n=0}
        {
            has_vw=0;
            for(i=12;i<=NF;i++){
                if($i ~ /^vW:i:/){has_vw=1}
            }
            if(has_vw==0){n++}
        }
        END{print n}')

    VW1_IN=$(samtools view -@ "${THREADS}" "${IN_BAM}" | \
        awk 'BEGIN{n=0}
        {
            for(i=12;i<=NF;i++){
                if($i=="vW:i:1"){n++}
            }
        }
        END{print n}')

    VW_OTHER_IN=$(samtools view -@ "${THREADS}" "${IN_BAM}" | \
        awk 'BEGIN{n=0}
        {
            has_other=0;
            for(i=12;i<=NF;i++){
                if($i ~ /^vW:i:/ && $i!="vW:i:1"){has_other=1}
            }
            if(has_other==1){n++}
        }
        END{print n}')

    echo "  Total input alignments: ${TOTAL_IN}"
    echo "  Input alignments with no vW tag: ${NO_VW_IN}"
    echo "  Input alignments with vW=1: ${VW1_IN}"
    echo "  Input alignments with other vW tag: ${VW_OTHER_IN}"

    # filter BAM
    echo "  Filtering BAM..."

    samtools view -@ "${THREADS}" -h "${IN_BAM}" | \
    awk 'BEGIN{OFS="\t"}
        /^@/ {
            print;
            next;
        }
        {
            has_vw=0;
            pass_vw=0;

            for(i=12;i<=NF;i++){
                if($i ~ /^vW:i:/){
                    has_vw=1;
                    if($i=="vW:i:1"){
                        pass_vw=1;
                    }
                }
            }

            if(has_vw==0 || pass_vw==1){
                print;
            }
        }' | \
    samtools view -@ "${THREADS}" -b -o "${TMP_BAM}" -

    # sort and index filtered BAM
    echo "  Sorting filtered BAM..."

    samtools sort \
        -@ "${THREADS}" \
        -m 3G \
        -o "${OUT_BAM}" \
        "${TMP_BAM}"

    rm -f "${TMP_BAM}"

    echo "  Indexing filtered BAM..."

    samtools index -@ "${THREADS}" "${OUT_BAM}"

    # generate stats
    echo "  Generating stats..."

    local TOTAL_OUT
    local NO_VW_OUT
    local VW1_OUT
    local VW_OTHER_OUT

    TOTAL_OUT=$(samtools view -@ "${THREADS}" -c "${OUT_BAM}")

    samtools flagstat -@ "${THREADS}" "${IN_BAM}" > "${STATS_DIR}/${SAMPLE}.input.flagstat.txt"
    samtools flagstat -@ "${THREADS}" "${OUT_BAM}" > "${STATS_DIR}/${SAMPLE}.wasp_filtered.flagstat.txt"
    samtools idxstats "${OUT_BAM}" > "${STATS_DIR}/${SAMPLE}.wasp_filtered.idxstats.txt"

    NO_VW_OUT=$(samtools view -@ "${THREADS}" "${OUT_BAM}" | \
        awk 'BEGIN{n=0}
        {
            has_vw=0;
            for(i=12;i<=NF;i++){
                if($i ~ /^vW:i:/){has_vw=1}
            }
            if(has_vw==0){n++}
        }
        END{print n}')

    VW1_OUT=$(samtools view -@ "${THREADS}" "${OUT_BAM}" | \
        awk 'BEGIN{n=0}
        {
            for(i=12;i<=NF;i++){
                if($i=="vW:i:1"){n++}
            }
        }
        END{print n}')

    VW_OTHER_OUT=$(samtools view -@ "${THREADS}" "${OUT_BAM}" | \
        awk 'BEGIN{n=0}
        {
            has_other=0;
            for(i=12;i<=NF;i++){
                if($i ~ /^vW:i:/ && $i!="vW:i:1"){has_other=1}
            }
            if(has_other==1){n++}
        }
        END{print n}')

    # Report
    {
        echo "Sample: ${SAMPLE}"
        echo "Input BAM: ${IN_BAM}"
        echo "Output BAM: ${OUT_BAM}"
        echo
        echo "Input alignment counts:"
        echo "  Total input alignments: ${TOTAL_IN}"
        echo "  No vW tag: ${NO_VW_IN}"
        echo "  vW:i:1: ${VW1_IN}"
        echo "  Other vW tag: ${VW_OTHER_IN}"
        echo
        echo "Output alignment counts:"
        echo "  Total output alignments: ${TOTAL_OUT}"
        echo "  No vW tag: ${NO_VW_OUT}"
        echo "  vW:i:1: ${VW1_OUT}"
        echo "  Other vW tag: ${VW_OTHER_OUT}"
        echo
        echo "Removed alignments:"
        echo "  $((TOTAL_IN - TOTAL_OUT))"
        echo
        echo "Output BAM index:"
        echo "  ${OUT_BAM}.bai"
    } > "${REPORT}"

    echo "  Completed sample: ${SAMPLE}"
    echo "  End time: $(date)"
	echo "-----------------------------"
	echo ""
}

# Process assigned samples
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