#!/bin/bash
#SBATCH --job-name=t3-5
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=t3-5_post_filter_multiqc_%j.out
#SBATCH --error=t3-5_post_filter_multiqc_%j.err

# ----------------------------------------
# Perform summarize QC results for WASP-filtered BAMs using MultiQC
# 1. Summarize flagstat metrics after WASP-filtering. 
# 2. Calculate WASP retention metrics. 
# 3. Summarize flagstat and Picard metrics using MultiQC. 
#
# Usage:
#   sbatch t3-5_post_filter_multiqc.sh
# ----------------------------------------

set -euo pipefail

# load modules
module load conda/3-23.3.1
conda activate rnaseq

# Paths
BASE_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
ALIGN_DIR="${BASE_DIR}/data/3_STAR_Alignment"
STATS_DIR="${ALIGN_DIR}/stats"

SAMPLE_LIST="${BASE_DIR}/data/sample_list.txt"

OUT_DIR="${ALIGN_DIR}/multiqc_outputs"
PICARD_DIR="${OUT_DIR}/picard"

MULTIQC_OUT="${OUT_DIR}/multiqc_report"
MULTIQC_FILENAME="wasp_picard_multiqc.html"

FLAGSTAT_SUMMARY="${MULTIQC_OUT}/wasp_filter_flagstat_summary.tsv"
WASP_RETENTION_SUMMARY="${MULTIQC_OUT}/wasp_filter_retention_summary.tsv"

mkdir -p "${MULTIQC_OUT}"

# script info 
SCRIPT="t3-5_post_filter_multiqc.sh"
echo "******************************"
echo "Script name: ${SCRIPT}"
echo "SLURM job ID: ${SLURM_JOB_ID}"
echo "Start date and time: $(date)"
echo "******************************"
echo ""

# sample list 
if [[ ! -f "${SAMPLE_LIST}" ]]; then
    echo "Error: Sample list not found: ${SAMPLE_LIST}" 
    exit 1
fi

mapfile -t SAMPLES < "${SAMPLE_LIST}"

# summarizing flagstat metrics 
echo "1. Generating flagstat summary..."
echo ""

extract_flagstat_metrics(){

    local FILE="$1"

    if [[ ! -f "${FILE}" ]]; then
        echo -e "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"
        return 0
    fi

    awk '
        BEGIN {
            total = "NA"
            secondary = "NA"
            supplementary = "NA"
            duplicates = "NA"
            mapped = "NA"
            paired = "NA"
            read1 = "NA"
            read2 = "NA"
            properly_paired = "NA"
            with_itself_and_mate_mapped = "NA"
            singletons = "NA"
        }

        /in total/ {
            total = $1
        }

        /secondary/ {
            secondary = $1
        }

        /supplementary/ {
            supplementary = $1
        }

        /duplicates/ && $0 !~ /primary duplicates/ {
            duplicates = $1
        }

        / mapped \(/ && $0 !~ /mate mapped/ {
            mapped = $1
        }

        /paired in sequencing/ {
            paired = $1
        }

        /read1/ {
            read1 = $1
        }

        /read2/ {
            read2 = $1
        }

        /properly paired/ {
            properly_paired = $1
        }

        /with itself and mate mapped/ {
            with_itself_and_mate_mapped = $1
        }

        /singletons/ {
            singletons = $1
        }

        END {
            print total "\t" \
                  secondary "\t" \
                  supplementary "\t" \
                  duplicates "\t" \
                  mapped "\t" \
                  paired "\t" \
                  read1 "\t" \
                  read2 "\t" \
                  properly_paired "\t" \
                  with_itself_and_mate_mapped "\t" \
                  singletons
        }
    ' "${FILE}"
}

# write summary
{
    echo -e "sample\tstage\tflagstat_file\ttotal\tsecondary\tsupplementary\tduplicates\tmapped\tpaired_in_sequencing\tread1\tread2\tproperly_paired\twith_itself_and_mate_mapped\tsingletons"
} > "${FLAGSTAT_SUMMARY}"

for SAMPLE in "${SAMPLES[@]}"; do

    INPUT_FLAGSTAT="${STATS_DIR}/${SAMPLE}.input.flagstat.txt"
    WASP_FLAGSTAT="${STATS_DIR}/${SAMPLE}.wasp_filtered.flagstat.txt"
    RG_FLAGSTAT="${STATS_DIR}/${SAMPLE}.wasp_filtered.rg.flagstat.txt"

    METRICS="$(extract_flagstat_metrics "${INPUT_FLAGSTAT}")"
    if [[ ! -f "${INPUT_FLAGSTAT}" ]]; then
        echo "Missing input flagstat for ${SAMPLE}: ${INPUT_FLAGSTAT}" 
    fi
    echo -e "${SAMPLE}\tinput\t${INPUT_FLAGSTAT}\t${METRICS}" >> "${FLAGSTAT_SUMMARY}"

    METRICS="$(extract_flagstat_metrics "${WASP_FLAGSTAT}")"
    if [[ ! -f "${WASP_FLAGSTAT}" ]]; then
        echo "Missing WASP-filtered flagstat for ${SAMPLE}: ${WASP_FLAGSTAT}" 
    fi
    echo -e "${SAMPLE}\twasp_filtered\t${WASP_FLAGSTAT}\t${METRICS}" >> "${FLAGSTAT_SUMMARY}"

    METRICS="$(extract_flagstat_metrics "${RG_FLAGSTAT}")"
    if [[ ! -f "${RG_FLAGSTAT}" ]]; then
        echo "Missing read-group flagstat for ${SAMPLE}: ${RG_FLAGSTAT}" 
    fi
    echo -e "${SAMPLE}\twasp_filtered_rg\t${RG_FLAGSTAT}\t${METRICS}" >> "${FLAGSTAT_SUMMARY}"

done

echo "flagstat summary written to: ${FLAGSTAT_SUMMARY}"
echo ""

# 2. Calculate WASP retention stats 
echo "2. Summarizing WASP retention..."

awk -F'\t' '
    BEGIN {
        OFS = "\t"
        print "sample",
              "input_total",
              "wasp_filtered_total",
              "removed_total",
              "retention_fraction",
              "input_mapped",
              "wasp_filtered_mapped",
              "mapped_retention_fraction"
    }

    NR == 1 {
        next
    }

    {
        sample = $1
        stage = $2
        total = $4
        mapped = $8

        if (stage == "input") {
            input_total[sample] = total
            input_mapped[sample] = mapped
        }

        if (stage == "wasp_filtered") {
            wasp_total[sample] = total
            wasp_mapped[sample] = mapped
        }

        samples[sample] = 1
    }

    END {
        for (sample in samples) {

            it = input_total[sample]
            wt = wasp_total[sample]
            im = input_mapped[sample]
            wm = wasp_mapped[sample]

            if (it == "" || it == "NA") {
                it = "NA"
            }

            if (wt == "" || wt == "NA") {
                wt = "NA"
            }

            if (im == "" || im == "NA") {
                im = "NA"
            }

            if (wm == "" || wm == "NA") {
                wm = "NA"
            }

            removed = "NA"
            retention = "NA"
            mapped_retention = "NA"

            if (it != "NA" && wt != "NA") {
                removed = it - wt
                if (it > 0) {
                    retention = wt / it
                }
            }

            if (im != "NA" && wm != "NA" && im > 0) {
                mapped_retention = wm / im
            }

            print sample,
                  it,
                  wt,
                  removed,
                  retention,
                  im,
                  wm,
                  mapped_retention
        }
    }
' "${FLAGSTAT_SUMMARY}" | sort -k1,1 > "${WASP_RETENTION_SUMMARY}"

echo "WASP retention summary written to: ${WASP_RETENTION_SUMMARY}"
echo ""

# 3. run MultiQC
echo "3. Running MultiQC..."

multiqc \
    "${PICARD_DIR}" \
    "${STATS_DIR}" \
    "${FLAGSTAT_SUMMARY}" \
    "${WASP_RETENTION_SUMMARY}" \
    --outdir "${MULTIQC_OUT}" \
    --filename "${MULTIQC_FILENAME}" \
    --force

echo "MultiQC completed."
echo "Report saved to: ${MULTIQC_OUT}/${MULTIQC_FILENAME}"
echo "End time: $(date)"