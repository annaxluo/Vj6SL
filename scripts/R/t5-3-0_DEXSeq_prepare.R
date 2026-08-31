# Analysis of exon bin usage
# 0. process DEXSeq outputs and sample information for statistical tests. 

library(tidyverse)
library(DEXSeq)
source("./utils.R")

# define paths and variables 
base_path <- getwd()
ver_str <- "20260202"

seed <- "1785449922"
set.seed(seed)

# parameters
# for DRIMSeq filtering 
n.small <- 4 # number of min replicates per condition for DRIMSeq filtering
min_gene_expr = 10
min_feature_expr = 10
min_feature_prop = 0.05

data_path <- file.path(base_path, "data", "5_DEXSeq")
flattend_gff <- file.path(base_path, "data", "5_DEXSeq", "annotations", 
                          "GRCm39.dexseq.gff")

output_path_base <- file.path(base_path, "outputs")
output_path <- file.path(output_path_base, paste0("outputs_t5-DEXSeq_", ver_str))
if(!dir.exists(output_path))
  dir.create(output_path)

# 1. process sample ids and conditions ----------------------------------
sample_ids <- c("sample1", "sample2", "sample3", "sample4") # mock 
genotype <- c("Control", "Het", "Control", "Het") # mock 
brain_region <- c("Cortex", "Cortex", "Cortex", "Cortex") # mock 

sample_info <- data.frame(sample_id = sample_ids, 
                          genotype = genotype, 
                          subject_id = sample_ids, 
                          brain_region = brain_region)

sample_info$genotype <- factor(sample_info$genotype, levels=c("Control", "Het"))
sample_info$subject_id <- factor(sample_info$subject_id, levels=sample_ids)

sample_info$brain_region <- factor(sample_info$brain_region, 
                                   levels=c("Cortex"))

# total reads 
mapped_reads <- read.delim(file.path(data_path, 
                                     "allele_expression_normalized", 
                                     "mapped_reads.tsv"))
sample_info <- merge(sample_info, mapped_reads, by="sample_id")

# save sample info 
out_fn1 <- file.path(output_path, "p1_sample_info.rds")
saveRDS(sample_info, out_fn1)

# 1. Read data, and apply DRIMSeq to filter lowly expressed features and genes----
count_files <- file.path(data_path, "5_DEXSeq", "counts", 
                         paste0(sample_info$sample_id, ".dexseq_counts.txt"))

names(count_files) <- sample_info$sample_id

counts_long <- bind_rows(lapply(names(count_files), function(ii){
  read_dexseq_count_file(count_files[[ii]], ii)
}))

# change to count matrix
counts_ <- counts_long %>%
  mutate(sample_id = as.character(sample_id),
         count = as.integer(count)) %>%
  tidyr::pivot_wider(id_cols = c(gene_id, feature_id),
                     names_from = sample_id,
                     values_from = count,
                     values_fill = 0) %>% 
  as.data.frame()

# prepare sample data for DRIMSeq 
samples_dm <- data.frame(sample_id = as.character(sample_info$sample_id), 
                         group = sample_info$genotype)

ddm_1 <- DRIMSeq::dmDSdata(counts = counts_,
                           samples = samples_dm)

# 2. Filter lowly expressed features ----------------------------------
ddm_2 <- DRIMSeq::dmFilter(ddm_1,
                           min_samps_gene_expr = n.small, 
                           min_samps_feature_expr = n.small, 
                           min_samps_feature_prop = n.small, 
                           min_gene_expr = min_gene_expr, 
                           min_feature_expr = min_feature_expr,
                           min_feature_prop = min_feature_prop)

# save object 
out_fn2 <- file.path(output_path, "p2_data-DRIMseq_filtering.rds")
saveRDS(ddm_2, out_fn2)

# 3. Create DEXSeq object for genotype (no random effect variable) --------------
library(org.Mm.eg.db) 

# make DEXSeq object -----------------
dxd_3 <- DEXSeq::DEXSeqDataSetFromHTSeq(
  countfiles = count_files,
  sampleData = sample_info,
  design = ~ sample + exon + brain_region:exon + genotype:exon,
  flattenedfile = flattend_gff
)

# add gene symbols 
gene_ids <-  lapply(as.character(rowData(dxd_3)$groupID), function(rr){
  unlist(strsplit(rr, "\\+"))
})

gene_ids_unique <- unique(do.call("c", gene_ids)) 

mapping_ <- AnnotationDbi::mapIds(
  org.Mm.eg.db,
  keys = gene_ids_unique,
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)

gene_symbols <- lapply(gene_ids, function(ll){
  paste(mapping_[ll], collapse = "+")
})

rowData(dxd_3)$gene_symbol <- do.call("c", gene_symbols)

# save object 
out_fn <- file.path(output_path, "p3_data-dex_3_genotype.rds")
saveRDS(dxd_3, out_fn)

# filter features and genes ----------------------------------
exon_bins_filtered <- DRIMSeq::counts(ddm_2)[,2]

# filter features 
feature_ids_dxd <- paste0(rowData(dxd_3)[,"groupID"], ":", 
                          sub("^E", "", rowData(dxd_3)[,"featureID"]))
keep_ <- feature_ids_dxd %in% exon_bins_filtered

dxd_3_filtered <- dxd_3[keep_,]

# remove monoexonic genes 
gene_bin_counts <- table(rowData(dxd_3_filtered)$groupID)
genes_keep <- names(gene_bin_counts)[gene_bin_counts >= 2]

dxd_3.1 <- dxd_3_filtered[rowData(dxd_3_filtered)$groupID %in% genes_keep, ]

out_fn3 <- file.path(output_path, "p3_data-dex_3-1_genotype.rds")
saveRDS(dxd_3.1, out_fn3)

# save the mapping 
out_fn4 <- file.path(output_path, "p4_gene_symbol_mapping.rds")
saveRDS(mapping_, out_fn4)

