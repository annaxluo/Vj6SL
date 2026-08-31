#!/usr/bin/env bash

# module setup for the vj6sl workflow.

module purge

# core bioinformatics tools
module load bcftools/1.18
module load fastqc/0.12.1
module load star/2.7.8a
module load gffread/0.12.7
module load htslib/1.18
module load samtools/1.18
module load gatk/4.5.0.0
module load subread/2.0.0
module load Salmon/1.10.1

module load conda/3-23.3.1

echo "Loaded modules:"
module list
