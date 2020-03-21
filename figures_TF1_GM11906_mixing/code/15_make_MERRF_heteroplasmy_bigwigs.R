library(dplyr)
library(GenomicRanges)
library(data.table)
library(rtracklayer)
library(SummarizedExperiment)

perc.rank <- function(x) trunc(rank(x))/length(x)

# Import mscATAc
SE <- readRDS("../../../mtscATACpaper_large_data_files/intermediate/GM11906_combinedSE.rds")

barcodes_all <- colnames(SE)
bc_low <- barcodes_all[colData(SE)$merrf_af < 0.1]
bc_mid <- barcodes_all[colData(SE)$merrf_af > 0.1 & colData(SE)$merrf_af < 0.6]
bc_hi <- barcodes_all[colData(SE)$merrf_af > 0.6]

frags <- rbind(fread("../../../mtscATACpaper_large_data_files/source/cellranger_output/Mix_Fix_1hr_fragments.tsv.gz"), 
               fread("../../../mtscATACpaper_large_data_files/source/cellranger_output/Mix_Fix_6hr_fragments.tsv.gz"))

# Setup and export low heteroplasmy
gr_low <- frags %>% setnames(c("chr", "start", "end", "bc", "pcrn")) %>% filter(bc %in% bc_low) %>%
  makeGRangesFromDataFrame()
reads_coverage <- coverage(gr_low)/length(gr_low)*1000000
export.bw(reads_coverage, con = paste0("../../../mtscATACpaper_large_data_files/intermediate/het_level_bigwig/Merrf_lowHet.bw"))

# Setup and export mid heteroplasmy
gr_mid <- frags %>% setnames(c("chr", "start", "end", "bc", "pcrn")) %>% filter(bc %in% bc_mid) %>%
  makeGRangesFromDataFrame()
reads_coverage <- coverage(gr_mid)/length(gr_mid)*1000000
export.bw(reads_coverage, con = paste0("../../../mtscATACpaper_large_data_files/intermediate/het_level_bigwig/Merrf_midHet.bw"))

# Setup and export high heteroplasmy
gr_high <- frags %>% setnames(c("chr", "start", "end", "bc", "pcrn")) %>% filter(bc %in% bc_hi) %>%
  makeGRangesFromDataFrame()
reads_coverage <- coverage(gr_high)/length(gr_high)*1000000
export.bw(reads_coverage, con = paste0("../../../mtscATACpaper_large_data_files/intermediate/het_level_bigwig/Merrf_highHet.bw"))


