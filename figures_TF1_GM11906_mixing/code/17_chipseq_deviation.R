library(diffloop)
library(chromVAR)
library(Matrix)
library(dplyr)
library(SummarizedExperiment)
library(BSgenome.Hsapiens.UCSC.hg19)
library(stringr)
library(BuenColors)
library(ggbeeswarm)
library(ggrastr)
library(BiocParallel)
register(MulticoreParam(2))

# Import files / make nice names
fileshort <- list.files("../data/gm12878_encode/", full.names = FALSE)
files <- list.files("../data/gm12878_encode/", full.names = TRUE)
names <- gsub("UniPk.narrowPeak.gz", "", gsub("wgEncodeAwgTfbs", "", fileshort)) %>%
  gsub(pattern = "Pcr1x", replacement = "" ) %>%
  gsub(pattern = "V0422111", replacement = "" ) %>%
  gsub(pattern = "Iggmu", replacement = "" ) %>%
  gsub(pattern = "Gm12878", replacement = "," ) %>%
  str_split_fixed(pattern = ",", 2) %>% data.frame() %>%
  pull(X2) %>% as.character()

se <- readRDS("../../../mtscATACpaper_large_data_files/intermediate/GM11906_combinedSE.rds")
peaks <- rowRanges(se)
se <- addGCBias(se, genome = BSgenome.Hsapiens.UCSC.hg19)

# Built match matrix
mat <- sapply(1:length(files), function(i){
  v <- as.numeric(1:length(peaks) %in% queryHits(findOverlaps(peaks, bedToGRanges(files[i]))))
  return(v)
})
colnames(mat) <- names
mmat <- Matrix(mat)
remove(mat)

chipseqdev <- computeDeviations(se, mmat)
motif_scores <- data.matrix(assays(chipseqdev)[["z"]])
het <- colData(se)$merrf_af
data.frame(cor = cor(t(motif_scores), het, method = "spearman", use = "pairwise.complete.obs"), transcription_factor = rownames(motif_scores)) %>%
  arrange(desc(cor)) %>% mutate(rank = 1:n()) -> odf

odf$color <- ifelse(odf$cor > 0.05, "up", ifelse(odf$cor < -0.3, "down", "ns"))
p1 <- ggplot(odf, aes(x = rank, y = cor, fill = color)) +
  geom_bar(stat = "identity", color = "black") + pretty_plot(fontsize = 8) + L_border() + 
  scale_fill_manual(values = c("firebrick","lightgrey", "dodgerblue3"))+
  theme(legend.position = "none") + labs(x = "Rank", y = "Spearman correlation")
cowplot::ggsave2(p1, file = "../plots/chipseq_correlation_plot.pdf", width = 2, height = 2)

plot_df <- data.frame(
  level = ifelse(het < 0.1, "alow", ifelse(het < 0.6, "bmid", "chigh")),
  CTCF = motif_scores["Ctcf",],
  MEF = motif_scores["Mef2csc13268V0416101",]
)

pB <- ggplot(plot_df, aes(x = level, y = MEF)) +
  geom_quasirandom(size = 0.2) +
  geom_boxplot(color = "firebrick", fill = NA, outlier.shape = NA) +
  pretty_plot(fontsize = 8) + L_border() + labs(x = "", y = "MEF2C Deviation Score")
cowplot::ggsave2(pB, file = "../plots/MEF2C_plot.pdf", width = 2, height = 2)
