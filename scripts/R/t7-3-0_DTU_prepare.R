# Analysis of transcript usage
# 0. process Salmon outputs and sample information for statistical tests. 

library(tidyverse)
library(DRIMSeq)
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

salmon_path <- file.path(base_path, "data", "7_Salmon", "quant")

output_path_base <- file.path(base_path, "outputs")
output_path <- file.path(output_path_base, paste0("outputs_t7-Salmon_", ver_str))
if(!dir.exists(output_path))
  dir.create(output_path)

# 1. make "transcript-to-gene" table for `txImport` -----------------------------
gtf_fn <- file.path(base_path, "ref_data", "refdata-gex-GRCm39-2024-A", 
                    "genes", "genes.gtf")

gtf <- rtracklayer::import(gtf_fn)

tx2gene <- as.data.frame(gtf) %>% 
  filter(type=="transcript") %>% 
  select(tx_id = transcript_id, 
         gene_id = gene_id, 
         gene_name = gene_name) %>% 
  distinct()

# save to tsv
out_fn <- file.path(output_path, "p0_tx2gene_mm39.tsv")
write.table(tx2gene, file = out_fn, sep = "\t", quote = FALSE, row.names = FALSE)

# 2. Read Salmon data into DRIMSeq objects --------------------------------

# read Salmon outputs
fn_list <- file.path(salmon_path, sample_info$sample_id, "quant.sf")
names(fn_list) <- sample_info$sample_id

# txi <- tximport::tximport(
#   fn_list,
#   type = "salmon",
#   txOut = TRUE,
#   countsFromAbundance = "no"
# )

txi <- tximport::tximport(
  fn_list,
  type = "salmon",
  txIn = TRUE, 
  txOut = TRUE,
  countsFromAbundance = "dtuScaledTPM", 
  tx2gene = tx2gene
)

# make DRIMSeq object 
counts_mat <- txi$counts
counts_mat <- counts_mat[rowSums(counts_mat) > 0,]

counts_df <- as.data.frame(counts_mat)
counts_df$feature_id <- rownames(counts_df)

d <- counts_df %>% 
  inner_join(tx2gene, by=c("feature_id"="tx_id")) %>% 
  select(gene_id, feature_id, all_of(sample_info$sample_id))

ddm_1 <- dmDSdata(counts = d, samples = sample_info)

out_fn2 <- file.path(output_path, "p2_data-ddm_1.rds")
saveRDS(ddm_1, out_fn2)

# filter 
ddm_2 <- dmFilter(
  ddm_1,
  min_samps_feature_expr = n.small,
  min_samps_feature_prop = n.small,
  min_samps_gene_expr = n.small,
  min_feature_expr = min_feature_expr,
  min_feature_prop = min_feature_prop,
  min_gene_expr = min_gene_expr
)

out_fn2.1 <- file.path(output_path, "p1_data-ddm_2.rds")
saveRDS(ddm_2, out_fn2.1)










