library(ggrastr)
library(BuenColors)

pbdf <- readRDS("../output/PBMCatac_SignacSeurat_labelTransfer.rds")

tb <- theme(legend.position = "none",
            panel.grid = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.background = element_blank()) 

colors = c("#1f77b4","#d62728","#2ca02c","#ff7f0e","#9467bd","#8c564b","#e377c2","#7f7f7f",
           "#bcbd22","#17becf","#ad494a","#e7ba52","#8ca252","#756bb1","#636363","#aec7e8",
           "#ff9896","#98df8a","#ffbb78","#c5b0d5","#c49c94","#f7b6d2","#c7c7c7","#dbdb8d",
           "#9edae5","#e7969c","#e7cb94","#c7e9c0","#de9ed6","#d9d9d9")

pu <- ggplot(pbdf, aes(x = UMAP_1, y = UMAP_2, color = ATAC_snn_res.0.5)) +
  geom_point(size = 1) + scale_color_manual(values = colors) +tb

ggsave(pu,
       filename = "../plots/raster_cluster_umap_pbmc.png",
       width = 8, height = 8, units = "in", dpi = 300)

