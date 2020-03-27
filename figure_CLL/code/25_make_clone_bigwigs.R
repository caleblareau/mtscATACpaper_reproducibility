library(data.table)
library(rtracklayer)
library(dplyr)
library(GenomicRanges)

"%ni%" <- Negate("%in%")

cdf <- readRDS("../output/PT1_clone_definition.rds")
frags <- fread("../../../mtscATACpaper_large_data_files/source/cellranger_output/CLL_PT1_CD19pos_v12-mtMask_fragments.tsv.gz", header = FALSE) %>% data.frame()

clusters <- as.character(unique(cdf$cluster_id))

sapply(clusters, function(cluster){
  
  possible_ids <- cdf %>% filter(cluster_id == cluster) %>% pull(cell_id) %>% as.character()
  cluster_gr <- frags %>% filter(V4 %in% possible_ids) %>%
    setnames(c("chr", "start", "end", "V4", "PCRn")) %>%
    makeGRangesFromDataFrame()
  
  reads_coverage <- coverage(cluster_gr)/length(cluster_gr)*1000000
  export.bw(reads_coverage, con = paste0("../../../mtscATACpaper_large_data_files/intermediate/cll_clone_bigwigs/CLL-Patient1-Bcell-c", as.character(cluster), ".bw"))
  cluster
}) -> bulk2

# From manual pruning...
clonal_bc <- cdf %>% filter(cluster_id %ni% c("01", "02", "05")) %>% pull(cell_id) %>% as.character()
cluster_gr1 <- frags %>% filter(V4 %in% clonal_bc) %>%
  setnames(c("chr", "start", "end", "V4", "PCRn")) %>%
  makeGRangesFromDataFrame()

reads_coverage <- coverage(cluster_gr1)/length(cluster_gr1)*1000000
export.bw(reads_coverage, con = paste0("../../../mtscATACpaper_large_data_files/intermediate/cll_clone_bigwigs/CLL-Patient1-Bcell-", as.character("other-clonal"), ".bw"))

