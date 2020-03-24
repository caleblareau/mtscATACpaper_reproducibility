library(GenomicRanges)
library(SummarizedExperiment)
library(data.table)
library(dplyr)
library(BuenColors)
library(Matrix)
perc.rank <- function(x) trunc(rank(x))/length(x)

set.seed(1)

#--
# Code adapted from https://github.com/buenrostrolab/dscATAC_analysis_code/blob/master/mousebrain/code/01_geneScores-revisions.R
#-- 

# Import gene bodies; restrict to TSS
gdf <- read.table("../data/hg19-tss.bed", stringsAsFactors = FALSE)
tss <- data.frame(chr = gdf$V1, gene = gdf$V4, stringsAsFactors = FALSE)
tss$tss <-  ifelse(gdf$V5 == "+", gdf$V3, gdf$V2)
tss$start <- ifelse(tss$tss - 50000 > 0, tss$tss - 50000, 0)
tss$stop <- tss$tss + 50000

tss_idx <- makeGRangesFromDataFrame(tss, keep.extra.columns = TRUE)

# Import ATAC peaks
SE <- readRDS("../../../mtscATACpaper_large_data_files/intermediate/GM11906_combinedSE.rds")
adf <- data.frame(rowRanges(SE))[,c(1,2,3)]; colnames(adf) <- c("chr", "start", "end")
adf$mp <- (adf$start + adf$end)/2
atacgranges <- makeGRangesFromDataFrame(adf, start.field = "mp", end.field = "mp")

#Overlap between ATAC peaks and Ranges linker
ov <- findOverlaps(atacgranges, tss_idx)

# Do distance decay for the weights
dist <- abs(mcols(tss_idx)$tss[subjectHits(ov)] - start(atacgranges)[queryHits(ov)])
exp_dist_model <- exp(-1*dist/5000)

# Prep an outcome matrix
m <- Matrix::sparseMatrix(i = c(queryHits(ov), length(atacgranges)),
                          j = c(subjectHits(ov), length(tss_idx)),
                          x = c(exp_dist_model,0))
colnames(m) <- gdf$V4
m <- m[,which(Matrix::colSums(m) != 0)]

# import counts
counts <- data.matrix(assays(readRDS("../../../mtscATACpaper_large_data_files/intermediate/GM11906_combinedSE.rds"))[["counts"]])
geneScores <- data.matrix(t(m) %*% counts)
geneScores <- t(t(geneScores)/colSums(geneScores)) * 10000

# Do a permutation test
het <- colData(SE)$merrf_af
data.frame(cor = cor(t(geneScores), het, method = "spearman"), gene = rownames(geneScores)) %>%
  arrange(desc(cor)) -> odf

perc.rank <- function(x) trunc(rank(x))/length(x)
plot_example_df <- data.frame(merrf = het, nr2f2 = geneScores["SENP5",], colData(SE))
plot_example_df$coverage_quantile <- perc.rank(plot_example_df$mean_cov)

p1 <- ggplot(plot_example_df, aes(x = merrf*100, y = nr2f2, col = coverage_quantile)) +
  geom_point(size = 0.4) + pretty_plot(fontsize = 8) + L_border() + 
  labs(y = "SENP5 gene score", x = "m.8344A>G Heteroplasmy (%)") +
  theme(legend.position = "none") + 
  scale_color_gradientn(colors = jdb_palette("solar_basic")) 

cowplot::ggsave(p1, file = "../plots/MERRFallele_SENP5.pdf", width = 1.7, height = 1.7)

# Visualize the global correlation patters
data.frame(cor = cor(t(geneScores), sample(het), method = "spearman"), gene = rownames(geneScores)) %>%
  arrange(desc(cor)) -> odfp

odf$rank <- 1:dim(odf)[1]
odfp$rank <- 1:dim(odfp)[1]

odf$color <- ifelse(odf$cor > 0.2, "up", ifelse(odf$cor < -0.2, "down", "ns"))
p1 <- ggplot(odf, aes(x = rank, y = cor, color = color)) +
  geom_point_rast(inherit.aes = FALSE, size = 0.4, aes(x = rank, y = cor), color = "grey", data = odfp, raster.dpi = 1000) +
  geom_point_rast(size = 0.4, raster.dpi = 1000) + pretty_plot(fontsize = 8) + L_border() + 
  geom_hline(yintercept = -0.2, linetype = 2) +
  scale_color_manual(values = c("firebrick","black", "dodgerblue3"))+
  geom_hline(yintercept = 0.2, linetype = 2) +
  theme(legend.position = "none") + labs(x = "Rank", y = "Spearman correlation")

cowplot::ggsave2(p1, file = "../plots/MERRFallele_geneScoreCor.pdf", width = 3, height = 1.7)

