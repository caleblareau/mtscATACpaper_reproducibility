library(dplyr)
library(data.table)
library(reshape2)

# Import double positive population
dp <- read.table("../output/double_positive_population.tsv", header = TRUE)
bc <- dp %>% pull(sample)

# Function to determine the number of counts per haplotype of the two main variants (8202 and 8344) 
import_haplotype_counts <- function(hour){
  
  dt <- fread(paste0("../data/Mix",as.character(hour),"h_hs.ase.tsv.gz"), col.names = c("chr", "pos", "base", "BQ", "barcode", "read"))
  
  # Find reads that contained both variants
  dt[dt$barcode %in% bc,c("read",  "pos")] %>%
    distinct() %>%
    group_by(read) %>%
    summarize(count = n()) %>% 
    filter(count == 2) %>% pull(read) -> reads_both
  
  # Determine the best nucleotide per variant based on base quality
  dt[read %in% reads_both, ] %>% 
    group_by(read, pos, barcode) %>%
    top_n(n = 1, wt = BQ) %>%
    dcast(read + barcode ~ pos, value.var = "base") %>%
    data.frame() %>% 
    group_by(X8202,X8344) %>%
    summarize(count = n()) -> count_df
  count_df
}

# Summarize the proportion of reads per haplotype
big_df <- rbind(
  import_haplotype_counts(1),
  import_haplotype_counts(6)
) %>% group_by(X8202,X8344) %>%
  summarize(count = sum(count))

big_df

# Now build stacked bar graph in illustrator with information above...
