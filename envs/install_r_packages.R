# Install CRAN and Bioconductor packages required by the R workflow scripts

cran_pkgs <- c(
  "tidyverse",
  "data.table",
  "rstatix",
  "coin"
)

install.packages(
  setdiff(cran_pkgs, rownames(installed.packages())),
  repos = "https://cloud.r-project.org"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

bioc_pkgs <- c(
  "DEXSeq",
  "DRIMSeq",
  "tximport",
  "rtracklayer",
  "org.Mm.eg.db",
  "AnnotationDbi",
  "BiocParallel"
)

BiocManager::install(
  setdiff(bioc_pkgs, rownames(installed.packages())),
  ask = FALSE,
  update = FALSE
)