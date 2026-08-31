# Allele-specific expression
# 0. process ASEReadCounter outputs and sample information for statistical tests. 

library(tidyverse) 

# set paths and variables -----------------------------------
base_path <- getwd()
ver_str <- "20260201"

seed <- "1785449922"
set.seed(seed)

# paths
data_path <- file.path(base_path, "data", "4_ASE_Analysis")

output_path_base <- file.path(base_path, "outputs")
output_path <- file.path(output_path_base, paste0("outputs_t4-ASE_", ver_str))
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

# 2. read counts -----------------------------
ase_df <- do.call("rbind", lapply(seq(nrow(sample_info)), function(ii){
  df_ <- read.delim(file.path(data_path, 
                              "allele_counts_Q20", 
                              paste0(sample_info$sample_id[[ii]], ".ASEReadCounter.tsv"))) %>% 
    add_column(sample_id = sample_info$sample_id[[ii]], .before = 1) 
  df_ 
}))

ase_df <- merge(sample_info, ase_df, by="sample_id")

# compute normalized reads 
ase_df <- ase_df %>% 
  mutate(refCount.CPM=refCount / mapped_reads * 1e6, 
         altCount.CPM=altCount / mapped_reads * 1e6, 
         totalCount.CPM=totalCount / mapped_reads * 1e6)

out_fn2 <- file.path(output_path, "p2_ASE_counts.rds")
saveRDS(ase_df, out_fn2)

