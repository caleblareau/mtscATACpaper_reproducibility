library(dplyr)
library(data.table)
library(reshape2)

dp <- read.table("../output/double_positive_population.tsv", header = TRUE)
bc <- dp %>% pull(sample)

import_haplotype_counts <- function(hour){
  
  dt <- fread(paste0("zcat < ../data/input/Mix",as.character(hour),"h_hs.ase.tsv.gz"), col.names = c("chr", "pos", "base", "BQ", "barcode", "read"))
  dt[dt$barcode %in% bc,c("read",  "pos")] %>%
    distinct() %>%
    group_by(read) %>%
    summarize(count = n()) %>% 
    filter(count == 2) %>% pull(read) -> reads_both
  
  dt[read %in% reads_both, ] %>% 
    group_by(read, pos, barcode) %>%
    top_n(n = 1, wt = BQ) %>%
    dcast(read + barcode ~ pos, value.var = "base") %>%
    data.frame() %>% 
    group_by(X8202,X8344) %>%
    summarize(count = n()) -> count_df
  count_df
}

big_df <- rbind(
  import_haplotype_counts(1),
  import_haplotype_counts(6)
) %>% group_by(X8202,X8344) %>%
  summarize(count = sum(count))

big_df
