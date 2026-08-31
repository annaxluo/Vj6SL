# Software environments

This directory contains files for installing Python and R dependencies.

## Python environment

The Python helper scripts can be run using a conda/mamba environment.

```bash
module load conda/3-23.3.1

mamba env create -f envs/rnaseq_python.yaml
# or:
# conda env create -f envs/rnaseq_python.yaml

conda activate rnaseq
```

The required Python packages include:

- pandas
- pysam
- htseq


## R packages

Install R dependencies with:

```bash
module load R/4.3.0
Rscript envs/install_r_packages.R
```

The required packages include:

- DEXSeq
- DRIMSeq
- tximport
- rtracklayer
- org.Mm.eg.db
- AnnotationDbi
- BiocParallel
- tidyverse
- data.table
- rstatix
- coin
