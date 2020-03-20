library(data.table)
library(SummarizedExperiment)
library(dplyr)
library(Matrix)

# Function that processes per-base coverage for top barcodes
process_coverage_long <- function(mgatk_file, chem ){
  SE <- readRDS(mgatk_file)
  depth_vec <- colData(SE)$depth
  cutoff <- head(sort(depth_vec, decreasing = TRUE), 500) %>% tail(1)
  rm <- rowMeans(assays(SE)[["coverage"]][,depth_vec >= cutoff])
  return(data.frame(
    chem = chem,
    pos = 1:length(rm),
    cov = rm
  ))
}

# Handle the mgatk rds files in order to parse the coverage
path <- "../../../mtscATACpaper_large_data_files/source/other/scRNAseq_mgatk_benchmark/"
odf <- rbind(
  process_coverage_long(paste0(path,"pbmc-v5p-8k-mgatk.rds"), "10x-5p"),
  process_coverage_long(paste0(path,"pbmc-v3p-10k-mgatk.rds"), "10x-3pv3"),
  process_coverage_long(paste0(path,"pbmc-v2-8k-mgatk.rds"), "10x-3pv2"),
  process_coverage_long(paste0(path,"smartseq2-colonies-mgatk.rds"), "SmartSeq2")
)

# Visualize the per bp coverage
p1 <- ggplot(odf, aes(x = pos, y = cov, color = chem)) +
  geom_line(size = 0.5) +
  pretty_plot(fontsize = 8) + L_border() +
  scale_color_manual(values = c("dodgerblue", "firebrick", "purple3", "black")) +
  labs(x = "Position in mitochondrial genome", y = "coverage", color = "Technology")
cowplot::ggsave2(p1, file = "../plots/coverage_compare_bp_paper.pdf", width = 3.3, height = 2)

