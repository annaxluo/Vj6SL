#!/bin/bash
#SBATCH --job-name=t4-0
#SBATCH --output=t4-0_ase_prepare_%j.out
#SBATCH --error=t4-0_ase_prepare_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

# ----------------------------------------
# Prepare VCF files for ASEReadCounter  
# 1. Compress VCF file and create index.
# 2. Create reference dictionary.
# 
# Usage:
#   sbatch t4-0_ase_prepare.sh
# ----------------------------------------

set -euo pipefail

set +u
module load gatk/4.5.0.0
module load htslib/1.18 
module load bcftools/1.18
module load samtools/1.18
set -u

# paths 
BASE_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
VCF_DIR="${BASE_DIR}/variant_files"

VCF_IN="${VCF_DIR}/variant_list_mm39.vcf"
VCF_SORTED="${VCF_DIR}/variant_list_mm39.sorted.vcf"
VCF_GZ="${VCF_DIR}/variant_list_mm39.sorted.vcf.gz"

REF_FASTA="${BASE_DIR}/ref_data/refdata-gex-GRCm39-2024-A/fasta/genome.fa"

# script info 
SCRIPT="t4-0_ase_prepare.sh"
echo "******************************"
echo "Script name: ${SCRIPT}"
echo "SLURM job ID: ${SLURM_JOB_ID:-local}"
echo "Start date and time: $(date)"
echo "******************************"
echo ""

# 1. Compress VCF file and create index. 

echo "1. Prepare VCF files..."

# sort vcf file
if [[ ! -f "${VCF_SORTED}" ]]; then
    echo "  Sorting VCF file..."
	
	gatk SortVcf \
        -I "${VCF_IN}" \
        -O "${VCF_SORTED}"
fi

# compress
if [[ ! -f "${VCF_GZ}" ]]; then
    echo "  Compressing sorted VCF file..."
	bgzip -c "${VCF_IN}" > "${VCF_GZ}"
fi

# index
if [[ ! -f "${VCF_GZ}.tbi" ]]; then
    echo "  Indexing sorted VCF file..."
	tabix -p vcf "${VCF_GZ}"
fi

echo "  Done. "
echo ""

# 2. Create reference dictionary.
echo "2. Create reference dictionary..."

if [[ ! -f "${REF_FASTA}" ]]; then
    echo "  Error: reference FASTA not found: ${REF_FASTA}" 
    exit 1
fi

if [[ ! -f "${REF_FASTA}.fai" ]]; then
    echo "  Reference FASTA index not found. Creating..."
    samtools faidx "${REF_FASTA}"
fi

REF_DICT="${REF_FASTA%.*}.dict"
if [[ ! -f "${REF_DICT}" ]]; then
    echo "  Create reference dictionary..."
    gatk CreateSequenceDictionary \
        -R "${REF_FASTA}" \
        -O "${REF_DICT}"
fi

echo "  Done."
echo "End time: $(date)"