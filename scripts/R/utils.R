# utilies for bulk RNA-seq analysis 

# read DEXSeq counts 
read_dexseq_count_file <- function(fn, sample_id){
  
  x <- read.table(fn, header = FALSE, sep = "\t", stringsAsFactors = FALSE, 
                  col.names = c("exon_bin", "count"))
  
  x <- x[!grepl("^_", x$exon_bin), ]
  
  x$gene_id <- sub(":[^:]+$", "", x$exon_bin)
  x$feature_id <- x$exon_bin
  x$sample_id <- sample_id
  
  x[, c("gene_id", "feature_id", "sample_id", "count")]
}