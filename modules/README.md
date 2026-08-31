# HPC module installation. 

This directory contains installation scripts and Lmod Lua modulefiles for
software used by the workflow.

Set `SOFTWARE_ROOT` before installation.

Example:

```bash
export SOFTWARE_ROOT=/path/to/shared/software
bash modules/samtools/install_samtools_1.18.sh
```

After installing the software, copy or symlink the corresponding Lua files into
your system's Lmod module tree.

Example module tree:

```text
/path/to/modulefiles/
├── samtools/1.18.lua
├── htslib/1.18.lua
├── bcftools/1.18.lua
├── fastqc/0.12.1.lua
├── star/2.7.8a.lua
├── trimmomatic/0.39.lua
├── gatk/4.5.0.0.lua
├── gffread/0.12.7.lua
├── subread/2.0.0.lua
└── salmon/1.10.1.lua
```

## Software

| Tool | Version | Purpose |
|---|---:|---|
| FastQC | 0.12.1 | FASTQ quality control |
| Trimmomatic | 0.39 | Adapter and quality trimming |
| STAR | 2.7.8a | RNA-seq alignment |
| Samtools | 1.18 | BAM/SAM processing |
| HTSlib | 1.18 | bgzip/tabix and HTS libraries |
| BCFtools | 1.18 | VCF/BCF processing |
| GATK | 4.5.0.0 | Variant-aware processing and ASEReadCounter |
| gffread | 0.12.7 | GTF/GFF processing and transcript extraction |
| Subread | 2.0.0 | featureCounts |
| Salmon | 1.10.1 | Transcript quantification |

