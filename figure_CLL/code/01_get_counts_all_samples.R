library(SummarizedExperiment)
library(data.table)
library(GenomicRanges)
library(ggrastr)
library(dplyr)
library(BuenColors)

peaks <- diffloop::bedToGRanges("CLL_heme.fixedwidthpeaks.bed")

# function to get counts
getCountsFromFrags <- function(frag_gz_file,
                               peaks_gr,
                               barcodes){
  
  # Make GRanges of fragments that are solid for the cells that we care about
  frags_valid <- data.table::fread(paste0("zcat < ", frag_gz_file)) %>% 
    data.frame() %>% filter(V4 %in% barcodes) %>%  # filter for barcodes in our search set
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
  colData <- data.frame(
    sample = barcodes_found,
    depth = as.numeric(denom),
    FRIP = Matrix::colSums(m)/as.numeric(denom)
  )
  # Make sure that the SE can be correctly constructed
  stopifnot(all(colData$sample == colnames(m)))
  
  # Make summarized Experiment
  SE <- SummarizedExperiment::SummarizedExperiment(
    rowRanges = peaks_gr,
    assays = list(counts = m),
    colData = colData
  )
  return(SE)
}


do_all <- function(barcodes_file, fragments_file, short){
  
  # Import PBMC stuff
  bc_pub_pbmc <- as.character(read.table(barcodes_file)[,1])
  pbmc_SE <- getCountsFromFrags(fragments_file, peaks, bc_pub_pbmc)
  
  df <- data.frame(colData(pbmc_SE))
  df$density <- get_density(log10(df$depth + 1), df$FRIP*100)
  p1 <- ggplot(df %>% arrange(density), aes(x = log10(depth + 1),  y = FRIP*100, color = density)) +
    geom_point_rast(size = 0.1) + scale_color_gradientn(colors = jdb_palette("brewer_celsius")) +
    pretty_plot() + L_border() + labs(x = "log10 # Fragments", y = "FRIP") +
    theme(legend.position = "bottom") +
    geom_vline(xintercept = 3, linetype = 2) +
    geom_hline(yintercept = 50, linetype = 2)
  
  # Save outs
  cowplot::ggsave(p1, file = paste0("qc_plots/",short,"_FRIP_QC.pdf"), width = 2.5, height = 3)
  pbmc_SE2 <- pbmc_SE[, log10(df$depth) > 3 & df$FRIP > 0.5 ]
  saveRDS(pbmc_SE2, file = paste0("rds_SE/",short,"_QC.rds"))
}
do_all("public_data/atac_pbmc_5k_nextgem_barcodes.tsv", "public_data/atac_pbmc_5k_nextgem_fragments.tsv.gz", "PBMC_5k_nextgem")

do_all("public_data/atac_v1_pbmc_10k_barcodes.tsv", "public_data/atac_v1_pbmc_10k_fragments.tsv.gz", "PBMC_10k")
do_all("public_data/atac_pbmc_10k_nextgem_barcodes.tsv", "public_data/atac_pbmc_10k_nextgem_fragments.tsv.gz", "PBMC_10k_nextgem")
do_all("CLL_15min_CR-mtMask/outs/filtered_peak_bc_matrix/barcodes.tsv", "CLL_15min_CR-mtMask/outs/fragments.tsv.gz", "Patient1_CLL_15min")
do_all("CLL_5min_CR-mtMask/outs/filtered_peak_bc_matrix/barcodes.tsv", "CLL_5min_CR-mtMask/outs/fragments.tsv.gz", "Patient1_CLL_5min")
do_all("CLL_CD19pos_CR-mtMask/outs/filtered_peak_bc_matrix/barcodes.tsv", "CLL_CD19pos_CR-mtMask/outs/fragments.tsv.gz", "Patient1_CLL_pos")
do_all("CLL_CD19neg_CR-mtMask/outs/filtered_peak_bc_matrix/barcodes.tsv", "CLL_CD19neg_CR-mtMask/outs/fragments.tsv.gz", "Patient1_CLL_neg")

do_all("Patient2_CLL_CD19_pos_CR-mtMask/outs/filtered_peak_bc_matrix/barcodes.tsv", "Patient2_CLL_CD19_pos_CR-mtMask/outs/fragments.tsv.gz", "Patient2_CLL_pos")
do_all("Patient2_CLL_CD19_neg_CR-mtMask/outs/filtered_peak_bc_matrix/barcodes.tsv", "Patient2_CLL_CD19_neg_CR-mtMask/outs/fragments.tsv.gz", "Patient2_CLL_neg")
