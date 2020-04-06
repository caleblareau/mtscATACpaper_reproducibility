library(ggrastr)
library(BuenColors)
library(dplyr)

full_df <- readRDS("../output/21March2020_signac_process.rds")

tb <- theme(legend.position = "none",
            panel.grid = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.background = element_blank()) 

#---
# Make plots for main text
#---

pEpcam <- ggplot(full_df, aes(x = UMAP_1, y = UMAP_2, color = EPCAM)) +
  geom_point_rast(raster.dpi = 500, size = 2) + tb + scale_color_gradientn(colors =  c('lightgrey', 'blue'))

pCD45 <- ggplot(full_df, aes(x = UMAP_1, y = UMAP_2, color = PTPRC)) +
  geom_point_rast(raster.dpi = 500, size = 2) + tb + scale_color_gradientn(colors =  c('lightgrey', 'blue'))

pTREM1 <- ggplot(full_df, aes(x = UMAP_1, y = UMAP_2, color = TREM1)) +
  geom_point_rast(raster.dpi = 500, size = 2) + tb + scale_color_gradientn(colors =  c('lightgrey', 'blue'))

pClusters <- ggplot(full_df, aes(x = UMAP_1, y = UMAP_2, color = seurat_clusters)) +
  geom_point_rast(raster.dpi = 500, size = 2) + tb + scale_color_manual(values = c("dodgerblue3", "purple4", "forestgreen", "orange3", "firebrick", "black"))

#cowplot::ggsave2(cowplot::plot_grid(pEpcam,pCD45,pTREM1,pClusters, nrow = 1),
#                filename = "../plots/first_row.png",
#                width = 5.6, height = 1.4, units = "in", dpi = 1000)


#---
# Make plots for supplement
#---

pS1 <- ggplot(full_df, aes(x = UMAP_1, y = UMAP_2, color = IL1RL1)) +
  geom_point_rast(raster.dpi = 500, size = 2) + tb + scale_color_gradientn(colors =  c('lightgrey', 'blue'))

pS2 <- ggplot(full_df, aes(x = UMAP_1, y = UMAP_2, color = GATA3)) +
  geom_point_rast(raster.dpi = 500, size = 2) + tb + scale_color_gradientn(colors =  c('lightgrey', 'blue'))

pS3 <- ggplot(full_df, aes(x = UMAP_1, y = UMAP_2, color = KIAA0226L)) +
  geom_point_rast(raster.dpi = 500, size = 2) + tb + scale_color_gradientn(colors =  c('lightgrey', 'blue'))

pS4 <- ggplot(full_df, aes(x = UMAP_1, y = UMAP_2, color = KIT)) +
  geom_point_rast(raster.dpi = 500, size = 2) + tb + scale_color_gradientn(colors =  c('lightgrey', 'blue'))

#cowplot::ggsave2(cowplot::plot_grid(pS1, pS2, pS4, nrow = 1),
#                 filename = "../plots/supplemental_row.png",
#                 width = 4.5, height = 1.5, units = "in", dpi = 1000)



#---
# Make mito plots for supplement
#---
pM1 <- ggplot(shuf(full_df),aes(x = UMAP_1, y = UMAP_2, color = X16147C.T*100)) +
  geom_point_rast(raster.dpi = 500, size = 2) + tb + scale_color_gradientn(colors = c("lightgrey", "firebrick")) 

pM2 <- ggplot(full_df %>% arrange(X1227G.A), aes(x = UMAP_1, y = UMAP_2, color = X1227G.A*100)) +
  geom_point_rast(raster.dpi = 500, size = 2) + tb + scale_color_gradientn(colors = c("lightgrey", "firebrick")) 

pM3 <- ggplot(full_df %>% arrange(X9728C.T), aes(x = UMAP_1, y = UMAP_2, color = X9728C.T*100)) +
  geom_point_rast(raster.dpi = 500, size = 2) + tb + scale_color_gradientn(colors = c("lightgrey", "firebrick")) 

pM4 <- ggplot(full_df %>% arrange(X824T.C), aes(x = UMAP_1, y = UMAP_2, color = X824T.C*100)) +
  geom_point_rast(raster.dpi = 500, size = 2) + tb + scale_color_gradientn(colors = c("lightgrey", "firebrick")) 

pM5 <- ggplot(full_df %>% arrange(X3244G.A), aes(x = UMAP_1, y = UMAP_2, color = X3244G.A*100)) +
  geom_point_rast(raster.dpi = 500, size = 2) + tb + scale_color_gradientn(colors = c("lightgrey", "firebrick")) 


cowplot::ggsave2(cowplot::plot_grid(pM1, pM2, pM3, pM4, pM5, nrow = 1),
                 filename = "../plots/supplemental_row_mito_mutations.png",
                 width = 7.5, height = 1.5, units = "in", dpi = 1000)


