library(data.table)
library(SummarizedExperiment)
library(GenomicRanges)
library(dplyr)
library(ggrastr)
library(BuenColors)

peaks_gr <- diffloop::bedToGRanges("../data/merrf_tf1.fixedwidthpeaks.bed")

getCountsFromFrags_GM11906 <- function(frag_gz_file, meta_file){
  
  # Meta data
  meta <- read.table(meta_file, header = TRUE, stringsAsFactors = FALSE) %>%
    dplyr::filter(assign == "GM11906")
  
  barcodes <- meta %>% pull(cell_id)
  
  # Make GRanges of fragments that are solid for the cells that we care about
  frags_valid <- data.table::fread(cmd = paste0("zcat < ", frag_gz_file)) %>% 
    data.frame() %>% dplyr::filter(V4 %in% barcodes) %>%  # filter for barcodes in our search set
    GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "V1", start.field = "V2", end.field = "V3", keep.extra.columns = TRUE)
  
  # Get a denominator, per cell
  denom <- table(GenomicRanges::mcols(frags_valid)$V4)
  barcodes_found <- names(denom)
  
  # Get the overlaps with peaks
  ovPEAK <- GenomicRanges::findOverlaps(peaks_gr, frags_valid)
  
  # Establish a numeric index for the barcodes for sparse matrix purposes
  id <- factor(as.character(GenomicRanges::mcols(frags_valid)$V4), levels = barcodes_found)
  
  # Make sparse matrix with counts with peaks by  unique barcode
  countdf <- data.frame(peaks = S4Vectors::queryHits(ovPEAK),
                        sample = as.numeric(id)[S4Vectors::subjectHits(ovPEAK)]) %>%
    dplyr::group_by(peaks,sample) %>% dplyr::summarise(count = n()) %>% data.matrix()
  
  m <- Matrix::sparseMatrix(i = c(countdf[,1], length(peaks_gr)),
                            j = c(countdf[,2], length(barcodes_found)),
                            x = c(countdf[,3],0))
  colnames(m) <- barcodes_found
  
  # Make a polished colData
  colData_basic <- data.frame(
    cell_id = barcodes_found,
    depth = as.numeric(denom),
    FRIP = Matrix::colSums(m)/as.numeric(denom)
  )
  
  colData <- merge(colData_basic, meta, by = "cell_id")
  
  # Make summarized Experiment
  SE <- SummarizedExperiment::SummarizedExperiment(
    rowRanges = peaks_gr,
    assays = list(counts = m[,as.character(colData$cell_id)]),
    colData = colData
  )
  return(SE)
}

# import fragments and create a summarized experiment
SE1 <- getCountsFromFrags_GM11906("../../../mtscATACpaper_large_data_files/source/cellranger_output/Mix_Fix_1hr_fragments.tsv.gz", "../output/data1_meta.tsv")
SE6 <- getCountsFromFrags_GM11906("../../../mtscATACpaper_large_data_files/source/cellranger_output/Mix_Fix_6hr_fragments.tsv.gz", "../output/data6_meta.tsv")
SE <- cbind(SE1, SE6)
dim(SE)

# Now visualize QC

df <- data.frame(colData(SE))

df$log10_mtDNA <- log10(df$mean_cov + 1)
df$density <- get_density(df$log10_mtDNA, df$FRIP)

p1 <- ggplot(df %>% arrange(density), aes(x = mean_cov, y = FRIP*100, color = density)) +
  scale_x_log10(breaks = c(10, 50, 100, 500)) + 
  geom_point_rast(size = 0.3, raster.dpi = 600) + pretty_plot(fontsize = 8) + L_border() +
  geom_hline(yintercept = 40, linetype = 2) +
  geom_vline(xintercept = 50, linetype = 2) +
  labs(x = "mtDNA Coverage", y = "% nuclear reads in peaks") +
  scale_color_gradientn(colors = jdb_palette("solar_basic")) +
  theme(legend.position = "none")
cowplot::ggsave(p1, file = "../plots/GM11906_merrf_peaks_mtDNA.pdf", width = 1.8, height = 1.8)

SE_filt <- SE[,df$mean_cov >= 50 & df$FRIP > 0.4]
SE_filt <- chromVAR::filterPeaks(SE_filt, min_fragments_per_peak = 5)
saveRDS(SE_filt, file = "../../../mtscATACpaper_large_data_files/intermediate/GM11906_combinedSE.rds")
