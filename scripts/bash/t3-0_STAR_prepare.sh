#!/bin/bash
#SBATCH --job-name=t3-0
#SBATCH --output=t3-0_star_prepare_%j.out
#SBATCH --error=t3-0_star_prepare_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=4:00:00

# ----------------------------------------
# Prepare files for STAR+WASP alignment 
# 1. verify genome fasta and GFT files
# 2. generate STAR index for genes
# 3. verify vcf files 
# 4. verify trimmed fastq directory and create batches for job array 
# 
# Usage:
#   sbatch t3-0_STAR_prepare.sh
# ----------------------------------------

set -euo pipefail

# load modules
module load star/2.7.8a

# paths and variables
BASE_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
REF_DIR="${BASE_DIR}/ref_data/refdata-gex-GRCm39-2024-A"
GENOME_FASTA="${REF_DIR}/fasta/genome.fa"
GTF_FILE="${REF_DIR}/genes/genes.gtf.gz"
STAR_INDEX_DIR="${REF_DIR}/star_2.7.8a"

SAMPLE_LIST="${BASE_DIR}/data/sample_list.txt"
FASTQ_DIR="${BASE_DIR}/data/2_TrimmedData"
VCF_FILE="${BASE_DIR}/variant_files/variant_list_mm39.vcf"
OUTPUT_DIR="${BASE_DIR}/data/3_STAR_Alignment"

mkdir -p "${STAR_INDEX_DIR}" "${OUTPUT_DIR}"

# script info 
SCRIPT="t3-0_STAR_prepare.sh"
echo "******************************"
echo "Script name: ${SCRIPT}"
echo "SLURM job ID: ${SLURM_JOB_ID:-local}"
echo "Start date and time: $(date)"
echo "******************************"
echo ""

# task parameters
THREADS="${SLURM_CPUS_PER_TASK:-8}"

# 1. verify genome fasta and GFT files 
echo "1. Generating STAR Index."
echo "Start time: $(date)"

# verify input files
# gene GTF
if [[ ! -f "${GTF_FILE}" ]]; then
    echo "Error: GTF file not found: ${GTF_FILE}"
    exit 1
fi

if [[ "${GTF_FILE}" == *.gz ]]; then
    GTF_UNZIPPED="${REF_DIR}/genes/genes.gtf"
    if [[ ! -f "${GTF_UNZIPPED}" ]]; then
        gunzip -k "${GTF_FILE}"
    fi
    GTF_FILE="${GTF_UNZIPPED}"
fi

# genome fasta
if [[ ! -f "${GENOME_FASTA}" ]]; then
    echo "Error: Genome FASTA not found: ${GENOME_FASTA}"
    exit 1
fi

echo "Genome FASTA: ${GENOME_FASTA}"
echo "GTF file: ${GTF_FILE}"
echo "End time: $(date)"
echo ""

# 2. generate STAR index 
echo "2. Generating STAR index."
echo "Start time: $(date)"

STAR --runThreadN ${THREADS} \
     --runMode genomeGenerate \
     --genomeDir ${STAR_INDEX_DIR} \
     --genomeFastaFiles ${GENOME_FASTA} \
     --sjdbGTFfile ${GTF_FILE} \
     --sjdbOverhang 100

if [[ -f "${STAR_INDEX_DIR}/SA" && -f "${STAR_INDEX_DIR}/Genome" ]]; then
    echo "STAR index generated at: ${STAR_INDEX_DIR}"
else
    echo "Error: STAR index not found at: ${STAR_INDEX_DIR}"
    exit 1
fi

echo "STAR index directory: ${STAR_INDEX_DIR}"
echo "End time: $(date)"
echo ""

# 3. verify vcf files 
echo "3. Verifying vcf files."

if [[ -f "${VCF_FILE}" ]]; then
    echo "vcf file found: ${VCF_FILE}"
    VARIANT_COUNT=$(grep -v "^#" "${VCF_FILE}" | wc -l)
    echo "Number of variants: ${VARIANT_COUNT}"
else
    echo "Error: vcf file not found: ${VCF_FILE}"
    exit 1
fi
echo ""

# 4. verify trimmed fastq files  
echo "4. Verifying trimmed fastq files and creating batches."

if [[ ! -d "${FASTQ_DIR}" ]]; then
    echo "Error: FASTQ directory not found: ${FASTQ_DIR}"
    exit 1
fi

# create batches 
SAMPLE_COUNT=$(wc -l < "${SAMPLE_LIST}")

# check that all samples have both paired files
echo "Verifying paired files..."
MISSING_FILES=0
while read SAMPLE; do
    R1="${FASTQ_DIR}/${SAMPLE}_1_paired.fq.gz"
    R2="${FASTQ_DIR}/${SAMPLE}_2_paired.fq.gz"
    
    if [[ ! -f "${R1}" ]]; then
        echo "Missing R1: ${R1}"
        MISSING_FILES=$((MISSING_FILES + 1))
    fi
    if [[ ! -f "${R2}" ]]; then
        echo "Missing R2: ${R2}"
        MISSING_FILES=$((MISSING_FILES + 1))
    fi
done < "${SAMPLE_LIST}"

if [[ ${MISSING_FILES} -gt 0 ]]; then
    echo "Error: ${MISSING_FILES} files missing"
    exit 1
fi
echo "All paired files verified"

echo ""
echo "Verification complete. "
echo "End time: $(date)"
