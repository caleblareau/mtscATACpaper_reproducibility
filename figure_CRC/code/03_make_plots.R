library(ggrastr)
library(BuenColors)

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

cowplot::ggsave2(cowplot::plot_grid(pEpcam,pCD45,pTREM1,pClusters, nrow = 1),
                filename = "../plots/first_row.png",
                width = 5.6, height = 1.4, units = "in", dpi = 1000)


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

cowplot::ggsave2(cowplot::plot_grid(pS1, pS2, pS3, pS4, nrow = 1),
                 filename = "../plots/supplemental_row.png",
                 width = 6, height = 1.5, units = "in", dpi = 1000)


ggplot(mito_df,aes(x = UMAP_1, y = UMAP_2, color = X16147C.T*100)) +
  geom_point() + pretty_plot() + L_border() + scale_color_gradientn(colors = c("grey", "red")) +
  theme(legend.position = "bottom")

ggplot(mito_df,aes(x = UMAP_1, y = UMAP_2, color = X3244G.A*100)) +
  geom_point() + pretty_plot() + L_border() + scale_color_gradientn(colors = jdb_palette("solar_extra")) +
  theme(legend.position = "bottom")

#---
# Make heatmap of variants
#---
colVec <- c("dodgerblue3", "purple4", "forestgreen", "orange3", "firebrick", "orange", "green2")
names(colVec) <- as.character(0:6)

ha1 <- HeatmapAnnotation(df = data.frame(clusters = as.character(mito_df$seurat_clusters)), 
                         col = list(clusters = colVec)
)

af = data.matrix(assays(mut_se)[["allele_frequency"]][as.character(filter_df$variant),])
af[af>0.1] <- 0.1

pdf(paste0("../plots/heatmap_mtDNAvars.pdf"), width=4.5, height=4)
hm <- Heatmap(af, 
              col=as.character(jdb_palette("solar_rojos",type="continuous")),
              show_row_names = TRUE, 
              cluster_columns = TRUE,
              cluster_rows = TRUE,
              row_names_gp = gpar(fontsize = 5),
              column_names_gp = gpar(fontsize = 5),
              top_annotation = ha1,
              show_column_names = FALSE)
hm
dev.off()


