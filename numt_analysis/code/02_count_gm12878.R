library(data.table)
library(dplyr)
library(GenomicRanges)
"%ni%" <- Negate("%in%")

# Import GM peaks
gr <- fread("../data/encode_roadmap_peaks/wgEncodeAwgDnaseUwdukeGm12878UniPk.narrowPeak.gz") %>%
  makeGRangesFromDataFrame(seqnames.field = "V1", start.field = "V2", end.field = "V3")

# Get GM cells from the 10x data
gm_barcodes <- fread("../data/scatac/atac_hgmm_500_v1_singlecell.csv") %>% filter(is_hg19_cell_barcode == 1) %>% pull(barcode)
filt_frag <- fread("../data/scatac/atac_hgmm_500_v1_fragments.tsv.gz") %>% data.frame() %>% filter(V4 %in% gm_barcodes) %>%
  mutate(chr = gsub("hg19_", "", V1)) %>%
  makeGRangesFromDataFrame(seqnames.field = "chr", start.field = "V2", end.field = "V3", keep.extra.columns = TRUE)
ov <- findOverlaps(filt_frag, gr)
bc_counts <- table(mcols(filt_frag)$V4[queryHits(ov)])
summary(as.numeric(bc_counts/length(gr))) # 4% of peaks have reads

# Same code as script 01 to subset to mito blacklist
bl <- data.frame(fread("../data/hg19-blacklist.v2.bed.gz", col.names = c("chr", "start", "end", "r"))) %>%
  makeGRangesFromDataFrame()
mito <- data.frame(fread("../data/mito_blacklist_hg19.bed")) %>%
  makeGRangesFromDataFrame(seqnames.field = "V1", start.field = "V2", end.field = "V3")
mito_ss <- mito[1:length(mito) %ni% subjectHits(findOverlaps(bl, mito))]

# Get GM cells
gm_barcodes <- fread("../data/scatac/atac_hgmm_500_v1_singlecell.csv") %>% filter(is_hg19_cell_barcode == 1) %>% pull(barcode)
filt_frag <- fread("../data/scatac/atac_hgmm_500_v1_fragments.tsv.gz") %>% data.frame() %>% filter(V4 %in% gm_barcodes) %>%
  mutate(chr = gsub("hg19_", "", V1)) %>%
  makeGRangesFromDataFrame(seqnames.field = "chr", start.field = "V2", end.field = "V3", keep.extra.columns = TRUE)
ov2 <- findOverlaps(filt_frag, mito_ss, type = "within")
bc_counts_NUMT <- table(mcols(filt_frag)$V4[queryHits(ov2)])
summary(as.numeric(bc_counts_NUMT))
