# Analysis of exon bin usage
# 1. perform statistical tests to determine differential exon usage between 
#    genotypes for each brain region. 

library(tidyverse)
library(DEXSeq)
source("./utils.R") 

# paths -----------------------------------------------
base_path <- getwd()
ver_str <- "20260202"

seed <- "1785449922"
set.seed(seed)

output_path_base <- file.path(base_path, "outputs")
output_path <- file.path(output_path_base, paste0("outputs_t5-DEXSeq_", ver_str))

flattend_gff <- file.path(base_path, "data", "5_DEXSeq", "annotations", 
                          "GRCm39.dexseq.gff")

# parameters 
# for DRIMSeq filtering 
n.small <- 4 # number of min replicates per condition for DRIMSeq filtering
min_gene_expr = 10
min_feature_expr = 10
min_feature_prop = 0.05

# system -----------------------------------
args <- commandArgs(trailingOnly = TRUE)

if(!length(args)==1 && args[1] %in% c("Control", "Het")){
  stop(paste0("invalid input: ", args[1]))
}
used_genotype <- args[1]

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
BPPARAM <- BiocParallel::MulticoreParam(workers=n_cores, progressbar=TRUE) 

# for all regions --------------------------------------
sample_table <- readRDS(file.path(output_path, "p0_sample_info.rds"))

# mapping
mapping_df <- readRDS(file.path(output_path, "p4_gene_symbol_mapping.rds"))

# process data for each brain region -----------------------------
res_output_path <- file.path(output_path, paste0("p5_DEU-genotype_", used_genotype))
if(!dir.exists(res_output_path)) 
  dir.create(res_output_path)

sample_table_sub <- sample_table %>% filter(genotype==used_genotype) %>% 
  droplevels()

count_files_sub <- file.path(data_path, "5_DEXSeq", "counts", 
                             paste0(sample_table_sub$sample_id, ".dexseq_counts.txt"))
names(count_files_sub) <- sample_table_sub$sample_id

# 1. read DRIMSeq object and filter features --------------------------
print("1. Preparing DRIMSeq object...\n")

counts_long <- bind_rows(lapply(names(count_files_sub), function(ii){
  read_dexseq_count_file(count_files_sub[[ii]], ii)
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
samples_dm <- data.frame(sample_id = as.character(sample_table_sub$sample_id), 
                         group = sample_table_sub$genotype)

ddm_1.0 <- DRIMSeq::dmDSdata(
  counts = counts_,
  samples = samples_dm
)

ddm_1.1 <- DRIMSeq::dmFilter(ddm_1.0,
                             min_samps_gene_expr = n.small, 
                             min_samps_feature_expr = n.small, 
                             min_samps_feature_prop = n.small, 
                             min_gene_expr = 10, 
                             min_feature_expr = 10,
                             min_feature_prop = 0.025)

# save object 
out_fn1.1 <- file.path(res_output_path, "p5-1_data-DRIMseq_filtering.rds")
saveRDS(ddm_1.1, out_fn1.1)

print(paste0("Saved DRIMSeq object to: ", out_fn1.1, "\n"))

# 2. make DEXSeq object -------------------------------
print("2. Makig DEXSeq object...\n")

dxd_1.2 <- DEXSeq::DEXSeqDataSetFromHTSeq(
  countfiles = count_files_sub,
  sampleData = sample_table_sub,
  design = ~ sample + exon + mouse_id:exon + brain_region:exon, 
  flattenedfile = flattend_gff
)

# add gene symbols 
gene_symbols <- mapping_df$gene_symbol[match(rowData(dxd_1.2)$groupID, 
                                             mapping_df$gene_id)]

rowData(dxd_1.2)$gene_symbol <- gene_symbols

# save object 
out_fn1.2 <- file.path(res_output_path, 
                       paste0("p5-1_data-dex_1-2_", used_genotype, ".rds"))
saveRDS(dxd_1.2, out_fn1.2)

print(paste0("Saved pre-filtered DEXSeq object to: ", out_fn1.2, "\n"))

# filter features and genes ----------------------------------
exon_bins_filtered <- DRIMSeq::counts(ddm_1.1)[,2]

# filter features 
feature_ids_dxd <- paste0(rowData(dxd_1.2)[,"groupID"], ":", 
                          sub("^E", "", rowData(dxd_1.2)[,"featureID"]))
keep_ <- feature_ids_dxd %in% exon_bins_filtered

dxd_1.2_filtered <- dxd_1.2[keep_,]

# remove monoexonic genes 
gene_bin_counts <- table(rowData(dxd_1.2_filtered)$groupID)
genes_keep <- names(gene_bin_counts)[gene_bin_counts >= 2]

dxd_1.3 <- dxd_1.2_filtered[rowData(dxd_1.2_filtered)$groupID %in% genes_keep, ]

out_fn1.3 <- file.path(res_output_path, 
                       paste0("p5-1_data-dex_1-3_", used_genotype, ".rds"))
saveRDS(dxd_1.3, out_fn1.3)

print(paste0("Saved post-filtered DEXSeq object to: ", out_fn1.3, "\n"))

# 3. Perform DEU test --------------------------------
print("Performing DEU test...\n")

# size factor and dispersions
dxd_2.0 <- estimateSizeFactors(dxd_1.3)
dxd_2.1 <- estimateDispersions(dxd_2.0, fitType="local", BPPARAM=BPPARAM)

# plot diagnostic plot 
out_fn2.0 <- file.path(res_output_path, "p5-2_diagnostic-dispersions.pdf")
pdf(out_fn2.0, width = 7, height = 7)
plotDispEsts(dxd_2.1)
dev.off()

out_fn2.1 <- file.path(res_output_path, paste0("p5-2_data-dex_2-2_", used_genotype, ".rds"))
saveRDS(dxd_2.1, out_fn2.1)

# Perform DE test for brain region -----------------------------------------
res <- testForDEU(
  dxd_2.1,
  reducedModel = ~ sample + exon + mouse_id:exon,
  BPPARAM = BPPARAM
)

# fold change 
res <- estimateExonFoldChanges(
  res,
  fitExpToVar = "brain_region", 
  BPPARAM = BPPARAM
)

dxr <- DEXSeqResults(res)
dxr_df <- as.data.frame(dxr)

out_fn2.2.1 <- file.path(res_output_path, paste0("p5-2_DEU-res_", used_genotype, ".rds"))
outs <- list("res"=res, 
             "dxr"=dxr)
saveRDS(outs, out_fn2.2.1)

out_fn2.2.2 <- file.path(res_output_path, paste0("p5-2_DEU-res_df_", used_genotype, ".rds"))
saveRDS(dxr_df, out_fn2.2.2)

print("Done.")
