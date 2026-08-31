# Directory structure

```text
bulk-rnaseq-ase-deu-dtu/
├── README.md
├── structure.txt
├── config/
│   └── directory_structure.md
├── data/
│   ├── 0_RawData/
│   ├── 1_FastQC/
│   │   └── MultiQC/
│   ├── 2_TrimmedData/
│   ├── 3_STAR_Alignment/
│   │   ├── aligned_bams/
│   │   ├── filtered_bams/
│   │   ├── multiqc_outputs/
│   │   ├── reports/
│   │   └── stats/
│   ├── 4_ASE_Analysis/
│   │   ├── allele_counts_Q20/
│   │   ├── normalized_expression/
│   │   ├── reports/
│   │   ├── stats/
│   │   └── variant_read_sequences/
│   ├── 5_DEXSeq/
│   │   ├── annotation/
│   │   ├── count_qc/
│   │   ├── counts/
│   │   └── logs/
│   ├── 6_featureCounts/
│   │   ├── exon_counts/
│   │   ├── gene_counts/
│   │   ├── logs/
│   │   └── multiqc/
│   ├── 7_Salmon/
│   │   ├── index/
│   │   ├── logs/
│   │   ├── multiqc/
│   │   ├── quant/
│   │   └── reference/
│   ├── sample_list.txt
│   └── samplesheet.csv
├── envs/
│   ├── install_r_packages.R
│   ├── README.md
│   └── rnaseq_python.yaml
├── logs/
├── modules/
│   ├── bcftools/
│   ├── fastqc/
│   ├── gatk/
│   ├── gffread/
│   ├── htslib/
│   ├── salmon/
│   ├── samtools/
│   ├── star/
│   ├── subread/
│   ├── trimmomatic/
│   └── README.md
├── ref_data/
│   ├── GRCm39/
│   └── refdata-gex-GRCm39-2024-A/
│       ├── fasta/
│       ├── genes/
│       └── star_2.7.8a/
├── scripts/
│   ├── bash/
│   │   ├── load_modules.sh
│   │   ├── t0_prepare.sh
│   │   ├── t1_FastQC.sh
│   │   ├── t2_Trimmomatic.sh
│   │   ├── t3-0_STAR_prepare.sh
│   │   ├── t3-1_STAR_WASP_align.sh
│   │   ├── t3-2_WASP_filter.sh
│   │   ├── t3-3_add_readGroups.sh
│   │   ├── t3-4_picard_qc.sh
│   │   ├── t3-5_post_filter_multiqc.sh
│   │   ├── t4-0_ase_prepare.sh
│   │   ├── t4-1_ase_readcounter.sh
│   │   ├── t4-2_ase_extract_variant_reads.sh
│   │   ├── t4-3_ase_normalized_expr.sh
│   │   ├── t5-0_dexseq_prepare.sh
│   │   ├── t5-1_dexseq_count.sh
│   │   ├── t5-2_dexseq_summary.sh
│   │   ├── t6-1_featurecounts.sh
│   │   ├── t6-2_featurecounts_multiqc.sh
│   │   ├── t7-0_salmon_prepare.sh
│   │   ├── t7-1_salmon_quant.sh
│   │   └── t7-2_salmon_multiqc.sh
│   ├── python/
│   │   ├── DEXSeq/
│   │   │   ├── dexseq_count.py
│   │   │   └── dexseq_prepare_annotation.py
│   │   ├── t4_extract_variant_reads.py
│   │   └── t4_normalize_expression.py
│   └── R/
│       ├── t4-4-0_ASE_prepare.R
│       ├── t4-4-1_ASE_het-allelic-imbalance.R
│       ├── t4-4-2_ASE_allele-expr.R
│       ├── t5-3-0_DEXSeq_prepare.R
│       ├── t5-3-1_dexseq_DEU_by-genotype.R
│       ├── t7-3-0_DTU_prepare.R
│       └── utils.R
└── variant_files/
    └── variant_list_mm39.vcf
```

## Main directories

| Directory | Description |
|---|---|
| `data/` | Input and output data directories for each workflow step. Real data not tracked. |
| `envs/` | Conda/R environment files and dependency installation scripts. |
| `logs/` | SLURM log files. Logs are not tracked. |
| `modules/` | Scripts to install the required software and create modules on HPC. |
| `ref_data/` | Reference genome, annotation, and indexes. Not tracked. |
| `scripts/bash/` | SLURM/bash scripts for each workflow step. |
| `scripts/python/` | Python helper scripts for ASE and DEXSeq processing. |
| `scripts/R/` | R scripts for ASE, DEXSeq, and DTU analyses. |
| `variant_files/` | Mock/example VCF file location. Real variant files are not tracked. |


