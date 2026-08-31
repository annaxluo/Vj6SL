# Allele-specific expression
# 1. Analyze allelic imbalance in Het samples (using allele counts) binomial 
# tests at the group and individual levels. 

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

# read data ---------------------------------------------
min_REF_ALT_total <- 10 # remove samples with low total read counts

read_counts <- readRDS(file.path(output_path, "p2_ASE_counts.rds"))

read_counts_het <- read_counts %>% 
  filter(genotype=="Het") %>% 
  mutate(ALT.fraction = altCount / (refCount + altCount)) %>% 
  filter(refCount + altCount >= min_REF_ALT_total) # all samples used

read_counts_het <- droplevels(read_counts_het)

# 3. Perform group binomial test for each brain region -------------------------

stats_df3 <- lapply(levels(read_counts_het$brain_region), function(region_){
  
  df_ <- read_counts_het %>% filter(brain_region==region_)
  df_ <- droplevels(df_)
  
  res <- binom.test(
    x = sum(df_$altCount),
    n = sum(df_$refCount + df_$altCount),
    p = 0.5
  )
  
  data.frame(
    test = res$method, 
    alt_count = res$statistic, 
    all_count = res$parameter, 
    alt_fraction_estimate = res$estimate,
    p_value = res$p.value
  ) 
})

stats_df3 <- do.call("rbind", stats_df3)
rownames(stats_df3) <- levels(read_counts_het$brain_region)

out_fn3 <- file.path(output_path, "p3_stats_allele-counts_group.csv")
write.csv(stats_df3, out_fn3)


# 4. Perform binomial test on individual samples-------------------------------

stats_df4 <- lapply(levels(read_counts_het$brain_region), function(region_){
  
  ret_ <- do.call("rbind", lapply(levels(read_counts_het$subject_id), function(s_){
    df_ <- read_counts_het %>% filter(brain_region==region_ & subject_id==s_)
    df_ <- droplevels(df_)
    
    res <- binom.test(
      x = sum(df_$refCount),
      n = sum(df_$altCount + df_$refCount),
      p = 0.5
    )
    
    data.frame(
      sample_id = df_$sample_id, 
      genotype = df_$genotype, 
      subject_id = df_$subject_id, 
      brain_region = df_$brain_region, 
      test = res$method, 
      alt_count = res$statistic, 
      all_count = res$parameter, 
      alt_fraction_estimate = res$estimate,
      p_value = res$p.value
    ) 
  }))
  
  rownames(ret_) <- ret_$subject_id
  ret_ 
})

stats_df4 <- do.call("rbind", stats_df4)

out_fn4 <- file.path(output_path, "p4_stats_allele-counts_sample.csv")
write.csv(stats_df4, out_fn4)

# 


