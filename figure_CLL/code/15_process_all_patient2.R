library(SummarizedExperiment)
library(BuenColors)
library(dplyr)
library(Matrix)
library(BuenColors)
library(data.table)
library(ggbeeswarm)
library(matrixStats)
library(stringr)

source("../../global_functions/get_allele_freq_mat.R")

# Import projections and raw mtDNA  mutations
bcsdf <- fread("../output/projection_dfs/CLL_PT2_CD19neg_projection.tsv", header = TRUE) %>% data.frame()
SE <- readRDS("../../../mtscATACpaper_large_data_files/source/mgatk_output/CLL_PT2_CD19neg_v12-mtMask_mgatk.rds")
cd <- data.frame(colData(SE))
af <- computeAFMutMatrix(SE)

# Process/annotate mutations
df <- data.frame(barcode = colnames(af), (data.matrix(t(af[c("12980G>A", "4853G>A"),])) >= 0.2), depth = cd$depth) %>% filter(depth > 20)
mdfff <- merge(bcsdf, df, by.x = "barcode", by.y = "barcode")
mdfff %>% group_by(cluster) %>%
  summarize(m12980 = sum(X12980G.A), 
            m4853 = sum( X4853G.A ),
            n_cells_total = n()
  ) %>% data.frame()

# Make a visualization -- some extra variables for eash of plotting of mutations that will only be in scRNAseq
mdfff$score <- ifelse(mdfff$X12980G.A, "12980G>A", ifelse(mdfff$X4853G.A, "4853G>A", "znone"))
p1a <- ggplot(mdfff %>% arrange(desc(score)), aes(x = umap1, y =umap2, color = score)) +
  geom_point_rast(size = 0.2, raster.dpi = 1000) + pretty_plot(fontsize = 7) + L_border() +
  scale_color_manual(values = c("znone"="lightgrey", "LEF1" = "dodgerblue4","HCST" = "orange2", "4853G>A" = "dodgerblue", "12980G>A"="firebrick")) + 
  theme(legend.position = "none") + labs(x = "UMAP1", y = "UMAP2")

cowplot::ggsave2(p1a, file = "../plots/UMAP_CD19neg_CLL2_projection.pdf", 
                width = 1.7, height = 1.7)

