library(dplyr)
library(GenomicRanges)
library(data.table)
library(rtracklayer)

perc.rank <- function(x) trunc(rank(x))/length(x)

# Import mscATAc
SE <- readRDS("../data/big_gi/GM11906_combinedSE.rds")

barcodes_all <- colnames(SE)
bc_low <- barcodes_all[colData(SE)$merrf_af < 0.1]
bc_mid <- barcodes_all[colData(SE)$merrf_af > 0.1 & colData(SE)$merrf_af < 0.6]
bc_hi <- barcodes_all[colData(SE)$merrf_af > 0.6]

frags <- rbind(fread(cmd = "zcat < ../data/big_gi/Mix_Fix_1hr_fragments.tsv.gz"), 
               fread(cmd = "zcat < ../data/big_gi/Mix_Fix_6hr_fragments.tsv.gz"))

# Setup up low heteroplasmy
gr_low <- frags %>% setnames(c("chr", "start", "end", "bc", "pcrn")) %>% filter(bc %in% bc_low) %>%
  makeGRangesFromDataFrame()
reads_coverage <- coverage(gr_low)/length(gr_low)*1000000
export.bw(reads_coverage, con = paste0("../output/het_level_bigwig/Merrf_lowHet.bw"))

# Setup up mid heteroplasmy
gr_mid <- frags %>% setnames(c("chr", "start", "end", "bc", "pcrn")) %>% filter(bc %in% bc_mid) %>%
  makeGRangesFromDataFrame()
reads_coverage <- coverage(gr_mid)/length(gr_mid)*1000000
export.bw(reads_coverage, con = paste0("../output/het_level_bigwig/Merrf_midHet.bw"))

# Setup up high heteroplasmy
gr_high <- frags %>% setnames(c("chr", "start", "end", "bc", "pcrn")) %>% filter(bc %in% bc_hi) %>%
  makeGRangesFromDataFrame()
reads_coverage <- coverage(gr_high)/length(gr_high)*1000000
export.bw(reads_coverage, con = paste0("../output/het_level_bigwig/Merrf_highHet.bw"))


