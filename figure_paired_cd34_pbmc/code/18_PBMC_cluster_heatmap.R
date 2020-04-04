library(BuenColors)
library(dplyr)
library(reshape2)

# Make a data frame showing how atac clusters correspond to scRNA-seq clusters
df <- readRDS("../output/PBMC_clone_DF.rds")

long_df <- df %>% group_by(ATAC_snn_res.0.5, predicted.id) %>%
  summarize(count = n()) %>% ungroup() %>%
  group_by(ATAC_snn_res.0.5) %>% mutate(prop = count / sum(count))

reshape2::dcast(long_df, ATAC_snn_res.0.5 ~ predicted.id, value.var = "prop", fill = 0)

long_df$predicted.id <- factor(as.character(long_df$predicted.id), rev(unique(as.character(long_df$predicted.id))))
pOut <- ggplot(long_df, aes(y = predicted.id, x = ATAC_snn_res.0.5, fill = prop)) + 
  geom_tile() + pretty_plot(fontsize = 6) + L_border() +
  scale_fill_gradientn(colors = jdb_palette("solar_rojos")) +
  scale_y_discrete(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +
  theme(legend.position = "none") + labs(x = "ATAC Cluster", y = "scRNA transfer label")
cowplot::ggsave2(pOut, file = "../plots/heatmap_cluster_annotation.pdf", width = 2.1, height = 1.4)
