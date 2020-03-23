library(GenomicRanges)
library(dplyr)
library(data.table)
library(reshape2)
options(scipen=999)
"%ni%" <- Negate("%in%")

# Function to create a GRanges object based on resolution and step size
format_one <- function(chr, size, res = 10000000, by_res = 2000000){
  ranges_df <- data.frame(
    chr = chr,
    start = seq(from =0, to = size, by = by_res),
    end = seq(from =0, to = size, by = by_res) + res
  ) %>% mutate(idx = 1:n())
  
  # Make sure that at least half of the range is available based on max size
  ranges_df[size - ranges_df$end > by_res/2*(-1),]
  
}

# Setup cnv blocks
sizes <- read.table("chrom_hg19.sizes", stringsAsFactors = FALSE)
lapply(1:dim(sizes)[1], function(i){
  format_one(sizes[i,1], sizes[i,2])
}) %>% rbindlist() %>% data.frame() -> gr_df
region_gr <- makeGRangesFromDataFrame(gr_df)
labels <- paste(gr_df$chr, as.character(gr_df$start), as.character(gr_df$end), as.character(gr_df$idx), sep="_")

# Import fragments
make_basic_mat <- function(frag_file, barcodes){
  
  # Import and subset frags
  frags <- fread(frag_file)
  frags <- frags[V4 %in% barcodes]
  
  frag_gr <- makeGRangesFromDataFrame(data.frame(frags) %>% setnames(c("chr", "start", "end", "bc", "pcr")), keep.extra.columns = TRUE)
  
  # Remove blacklist
  bl <- diffloop::bedToGRanges("hg19.full.blacklist.bed")
  ov_bl <- findOverlaps(frag_gr, bl)
  boo <- (1:length(frag_gr) %ni% queryHits(ov_bl)) & width(frag_gr) < 2000
  print(sum(!boo))
  frag_gr <- frag_gr[boo]
  
  # Find overlaps
  ov <- findOverlaps(frag_gr, region_gr)
  barcodes <- mcols(frag_gr)[["bc"]]
  sum_df <- data.frame(barcode = barcodes[queryHits(ov)],
                       region_idx = subjectHits(ov)) %>% group_by(barcode, region_idx) %>% summarize(count = n())
  square <- reshape2::dcast(sum_df, region_idx ~ barcode, value.var = "count", fill = 0)
  row_labels <- labels[square[,1]]
  square <- square[,-1]
  rownames(square) <- row_labels
  
  # Define consistent set of regions
  regions_pbmc10k <- readRDS("vector_of_regions.rds")
  square <- square[regions_pbmc10k, ]
  
  square
}

