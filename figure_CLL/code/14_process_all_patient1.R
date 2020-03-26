library(SummarizedExperiment)
library(BuenColors)
library(dplyr)
library(Matrix)
library(BuenColors)
library(data.table)
library(ggbeeswarm)
library(dplyr)
library(matrixStats)
library(ComplexHeatmap)
library(circlize)
library(stringr)



computeAFMutMatrix <- function(SE){
  cov <- assays(SE)[["coverage"]]+ 0.001
  ref_allele <- as.character(rowRanges(SE)$refAllele)
  
  getMutMatrix <- function(letter){
    mat <- (assays(SE)[[paste0(letter, "_counts_fw")]] + assays(SE)[[paste0(letter, "_counts_rev")]]) / cov
    rownames(mat) <- paste0(as.character(1:dim(mat)[1]), "_", ref_allele, ">", letter)
    return(mat[ref_allele != letter & ref_allele != "N",])
  }
  
  rbind(getMutMatrix("A"), getMutMatrix("C"), getMutMatrix("G"), getMutMatrix("T"))
  
}


import_munge <- function(mgatk_file, proj_file){
  se <- readRDS(mgatk_file)
  af <- computeAFMutMatrix(se)
  
  bcsdf <- fread(proj_file, header = TRUE) %>% data.frame()
  df2 <- data.frame(barcode = colnames(af), data.matrix(t(af[c("5140_G>A", "1260_A>G", "14858_G>A", "1872_T>C"),] >= 0.2)), depth = colData(se)$depth) %>% filter(depth > 20)
  mdfff <- merge(bcsdf, df2, by.x = "barcode", by.y = "barcode")
  return(mdfff)
}

cd19neg <- import_munge("CLL_CD19neg_CR-mtMask_mgatk/final/CLL_CD19neg_CR-mtMask_mgatk.rds", "projection_dfs/Pt1_CD19neg_final.tsv")
min5 <- import_munge("CLL_15min_CR-mtMask_mgatk/final/CLL_15min_CR-mtMask_mgatk.rds", "projection_dfs/Pt1_15min_final.tsv")
min15 <- import_munge("CLL_5min_CR-mtMask_mgatk/final/CLL_5min_CR-mtMask_mgatk.rds", "projection_dfs/Pt1_5min_final.tsv")

df <- rbind(cd19neg, min5, min15)
df$score <- ifelse(df$X5140_G.A, "5140G>A", ifelse(df$X1260_A.G, "1260A>G", ifelse(df$X14858_G.A, "14858G>A", ifelse(df$X1872_T.C, "1872T>C", ""))))

p1a <- ggplot(rbind(df[df$score == "",], shuf(df[df$score != "",])), aes(x = umap1, y =umap2, color = score)) +
  geom_point_rast(size = 0.2, raster.dpi = 1000) + pretty_plot(fontsize = 7) + L_border() +
  scale_color_manual(values = c("lightgrey", "dodgerblue4", "dodgerblue", "orange2","firebrick")) + 
  theme(legend.position = "none") + labs(x = "UMAP1", y = "UMAP2")
cowplot::ggsave(p1a, file = "plots/UMAP_CD19neg_CLL1_projection.pdf", 
                width = 1.7, height = 1.7)

