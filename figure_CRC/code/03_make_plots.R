

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
#af[af>0.2] <- 0.2

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


