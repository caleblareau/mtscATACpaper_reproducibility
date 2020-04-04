library(data.table)
library(dplyr)
library(SummarizedExperiment)
library(Matrix)
library(BuenColors)
"%ni%" <- Negate("%in%")

cc = c("#1f77b4","#e377c2","#d62728","#ff7f0e","#9467bd","#8c564b","#7f7f7f",
       "#bcbd22","#17becf","#ad494a","#e7ba52","#2ca02c")

mgatk_se_cd34 <- readRDS("../output/filteredCD34_mgatk_calls.rds")
load("../output/CD34_umap_embedding_granja.rda")

tb <- theme(legend.position = "none",
            panel.grid = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.background = element_blank()) 

# Make plots of initial data
p1 <- ggplot(plot_df, aes(x= X1, y = X2, color = Clusters)) +
  geom_point_rast(size = 1, raster.dpi = 500) +
  labs(x = "UMAP1", y= "UMAP2", color = "") + tb +
  scale_color_manual(values = cc)

projection_df <- projection_df %>% filter(celltype != "Monocytes")
p2 <- ggplot(projection_df[dim(projection_df)[1]:1,], aes(x= umap1, y = umap2, color = celltype, label = celltype)) +
  geom_point_rast(size = 1, raster.dpi = 500) +
  labs(x = "UMAP1", y= "UMAP2", color = "C1 FACS ") +
  tb + 
  scale_color_manual(values = c(ejc_color_maps, "none" = "lightgrey", "Monocytes" = "orange2"))
cowplot::ggsave2(cowplot::plot_grid(p1, p2, nrow = 1),
                 filename = "../plots/raster_cluster_projection.png",
                 width = 2.8, height = 1.4, units = "in", dpi = 1000)

# visualize all variants
dm <- data.matrix(t(assays(mgatk_se_cd34)[["allele_frequency"]]))
dm[dm < 0.02] <- 0
big_plot_df <- data.frame(
  plot_df,
  dm
)


var_names <- colnames(data.frame(data.matrix(t(assays(mgatk_se_cd34)[["allele_frequency"]]))))

for(variant in var_names){
  print(variant)
  png(paste0("../plots/t02/",variant,"_02_threshold.png"), width = 4, height = 4.5, res = 300, units = "in")
  pp <- ggplot(big_plot_df %>% arrange((!!sym(variant))), aes_string(x = "X1", y = "X2", color = variant)) +
    geom_point() + scale_color_gradientn(colors = c("lightgrey", "red", "firebrick")) +
    pretty_plot() + L_border() + theme(legend.position = "bottom") + labs(x = "UMAP1", y = "UMAP2")
  print(pp)
  dev.off()
  
}




