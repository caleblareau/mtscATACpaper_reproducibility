library(BuenColors)
library(dplyr)
library(ggrastr)

# Import processed data
cd34_clone_df <- readRDS("../output/CD34_clone_DF.rds")
pbmc_clone_df <- readRDS("../output/PBMC_clone_DF.rds")
cd34_mut_se <- readRDS("../output/filteredCD34_mgatk_calls.rds")
pbmc_mut_se <- readRDS("../output/filteredpbmcs_mgatk_calls.rds")

tb <- theme(legend.position = "none",
            panel.grid = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.background = element_blank()) 


make_4plot_grid <- function(clone_one, variant, variant_name){
  
  # Make plots of clones
  cd34_clone_df$color_clone <- cd34_clone_df$mito_cluster == clone_one
  p_CD34_clone <- ggplot(cd34_clone_df %>% arrange(color_clone), aes(x= X1, y = X2, color = color_clone)) +
    geom_point_rast(size = 3, raster.dpi = 500) +
    tb + scale_color_manual(values = c("lightgrey", "dodgerblue3"))
  
  pbmc_clone_df$color_clone <- pbmc_clone_df$mito_cluster == clone_one
  p_PBMC_clone <- ggplot(pbmc_clone_df %>% arrange(color_clone), aes(x= UMAP_1, y = UMAP_2, color = color_clone)) +
    geom_point_rast(size = 3, raster.dpi = 500) +
    tb + scale_color_manual(values = c("lightgrey", "dodgerblue3"))
  
  cd34_clone_df$color_AF <- assays(cd34_mut_se)[["allele_frequency"]][variant,]
  p_CD34_AF <- ggplot(cd34_clone_df %>% arrange(color_AF), aes(x= X1, y = X2, color = color_AF)) +
    geom_point_rast(size = 3, raster.dpi = 500) +
    tb + scale_color_gradientn(colors = c("lightgrey", "firebrick"))
  
  pbmc_clone_df$color_AF <- assays(pbmc_mut_se)[["allele_frequency"]][variant,]
  p_PBMC_AF <- ggplot(pbmc_clone_df %>% arrange(color_AF), aes(x= UMAP_1, y = UMAP_2, color = color_AF)) +
    geom_point_rast(size = 3, raster.dpi = 500) +
    tb + scale_color_gradientn(colors = c("lightgrey", "firebrick"))
  
  cowplot::ggsave2(cowplot::plot_grid(p_CD34_AF,p_PBMC_AF,p_CD34_clone,p_PBMC_clone, nrow = 2),
                   filename = paste0("../plots/raster_clones_",variant_name,".png"),
                   width = 2.0, height = 2.0, units = "in", dpi = 500)
}

make_4plot_grid("119",  "12868G>A", "12868G-A")
make_4plot_grid("008",  "2788C>A", "2788C-A")
make_4plot_grid("032",  "3209A>G", "3209A-G")


sum(cd34_clone_df$mito_cluster == "119")
sum(cd34_clone_df$mito_cluster == "008")
sum(cd34_clone_df$mito_cluster == "032")

sum(pbmc_clone_df$mito_cluster == "119")
sum(pbmc_clone_df$mito_cluster == "008")
sum(pbmc_clone_df$mito_cluster == "032")
