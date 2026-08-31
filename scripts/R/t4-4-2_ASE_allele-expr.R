# Allele-specific expression
# 1. Analyze normalized allele expression levels using Wilcoxon tests. 

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

# read data ---------------------------------
min_REF_ALT_total <- 10 # remove samples with low total read counts

read_counts <- readRDS(file.path(output_path, "p2_ASE_counts.rds"))

read_counts_valid <- read_counts %>% 
  mutate(ALT.fraction = altCount / (refCount + altCount)) %>% 
  filter(refCount + altCount >= min_REF_ALT_total) # all samples used

read_counts_valid <- droplevels(read_counts_valid)


# 5. Expression of total transcripts: Control vs Het ------------------------
val_str <- "TOTAL_CPM"

stats_df5 <- lapply(levels(read_counts_valid$brain_region), function(region_){
  
  df_ <- read_counts_valid %>% filter(brain_region==region_)
  df_ <- droplevels(df_)
  
  colnames(df_)[colnames(df_)==val_str] <- "val"
  
  dat_test <- rstatix::wilcox_test(val ~ genotype, data = df_, 
                                   exact = TRUE)
  Z <- coin::wilcox_test(val  ~ genotype, data = df_)
  
  # effect size 
  effectsize <- rstatix::wilcox_effsize(val  ~ genotype, data = df_)
  
  group_means <- df_ %>% group_by(genotype) %>% 
    summarize(mean_val = mean(val))
  
  dat_test <- dat_test %>% data.table::as.data.table() %>% 
    mutate(Z = Z@statistic@teststatistic, 
           effsize = effectsize$effsize, 
           magnitude = effectsize$magnitude, 
           group1_mean = group_means$mean_val[[match(dat_test$group1,
                                                     group_means$genotype)]], 
           group2_mean = group_means$mean_val[[match(dat_test$group2,
                                                     group_means$genotype)]]) 
  
  dat_test <- dat_test %>% 
    add_column(brain_region = region_, .after = 1)
  
  dat_test[1,1] <- val_str
  dat_test
  
})

stats_df5 <- do.call("rbind", stats_df5)

out_fn5 <- file.path(output_path, "p5_stats_total-counts_by-genotype.csv")
write.csv(stats_df5, out_fn5)

# 6. Expression of REF allele: Control vs Het ---------------------------

stats_df6 <- lapply(levels(read_counts_valid$brain_region), function(region_){
  
  df_ <- read_counts_valid %>% filter(brain_region==region_)
  df_ <- droplevels(df_)
  
  dat_test <- rstatix::wilcox_test(REF_CPM ~ genotype, data = df_, 
                                   exact = TRUE)
  Z <- coin::wilcox_test(REF_CPM  ~ genotype, data = df_)
  
  # effect size 
  effectsize <- rstatix::wilcox_effsize(REF_CPM  ~ genotype, data = df_)
  
  group_means <- df_ %>% group_by(genotype) %>% 
    summarize(mean_val = mean(REF_CPM))
  
  dat_test <- dat_test %>% data.table::as.data.table() %>% 
    mutate(Z = Z@statistic@teststatistic, 
           effsize = effectsize$effsize, 
           magnitude = effectsize$magnitude, 
           group1_mean = group_means$mean_val[[match(dat_test$group1,
                                                     group_means$genotype)]], 
           group2_mean = group_means$mean_val[[match(dat_test$group2,
                                                     group_means$genotype)]]) 
  
  dat_test <- dat_test %>% 
    add_column(brain_region = region_, .after = 1)
  
  dat_test
  
})

stats_df6 <- do.call("rbind", stats_df6)

out_fn6 <- file.path(output_path, "p6_stats_REF-counts_by-genotype.csv")
write.csv(stats_df6, out_fn6)


# 7. 






