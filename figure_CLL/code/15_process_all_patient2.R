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

bcsdf <- fread("projection_dfs/Pt2_CD19neg_final.tsv", header = TRUE) %>% data.frame()

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
SE <- readRDS("Patient2_CLL_CD19_neg_CR-mtMask_mgatk/final/Patient2_CLL_CD19_neg_CR-mtMask_mgatk.rds")
cd <- data.frame(colData(SE))
af <- computeAFMutMatrix(SE)

# Import mutations
df <- data.frame(barcode = colnames(af), (data.matrix(t(af[c("12980_G>A", "4853_G>A"),])) >= 0.2), depth = cd$depth) %>% filter(depth > 20)
mdfff <- merge(bcsdf, df, by.x = "barcode", by.y = "barcode")
mdfff %>% group_by(cluster) %>%
  summarize(m12980 = sum(X12980_G.A), 
            m4853 = sum( X4853_G.A ),
            n_cells_total = n()
  ) %>% data.frame()

mdfff %>% filter(cluster != "Bcell" & (X12980_G.A | X4853_G.A))

mdfff$score <- ifelse(mdfff$X12980_G.A, "12980G>A", ifelse(mdfff$X4853_G.A, "4853G>A", "znone"))


p1a <- ggplot(mdfff %>% arrange(desc(score)), aes(x = umap1, y =umap2, color = score)) +
  geom_point_rast(size = 0.2, raster.dpi = 1000) + pretty_plot(fontsize = 7) + L_border() +
  scale_color_manual(values = c("znone"="lightgrey", "LEF1" = "dodgerblue4","HCST" = "orange2", "4853G>A" = "dodgerblue", "12980G>A"="firebrick")) + 
  theme(legend.position = "none") + labs(x = "UMAP1", y = "UMAP2")

cowplot::ggsave(p1a, file = "plots/UMAP_CD19neg_CLL2_projection.pdf", 
                width = 1.7, height = 1.7)

