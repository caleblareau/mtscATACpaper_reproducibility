library(data.table)
library(dplyr)
library(GenomicRanges)
"%ni%" <- Negate("%in%")

# Import blacklist of all regions and specifically mtDNA
bl <- data.frame(fread("../data/hg19-blacklist.v2.bed.gz", col.names = c("chr", "start", "end", "r"))) %>%
  makeGRangesFromDataFrame()
mito <- data.frame(fread("../data/mito_blacklist_hg19.bed")) %>%
  makeGRangesFromDataFrame(seqnames.field = "V1", start.field = "V2", end.field = "V3")
mito_ss <- mito[1:length(mito) %ni% subjectHits(findOverlaps(bl, mito))] # this is the mito-specific blacklist regions

# Import omnibus peak set
count_vec <- sapply(list.files("../data/encode_roadmap_peaks/", full.names = TRUE), function(file){
  gr <- fread(file) %>%
    makeGRangesFromDataFrame(seqnames.field = "V1", start.field = "V2", end.field = "V3")
  length(queryHits(findOverlaps(gr, mito_ss)))
})

# Estimate number of peaks affected per encode/roadmap peak set
summary(count_vec)
