library(data.table)
library(dplyr)
library(SummarizedExperiment)
library(BuenColors)

getNN <- function(mat_af, k.param = 10){
  set.seed(1)
  rownames(mat_af) <- make.unique(rownames(mat_af))
  obj <- FindNeighbors(mat_af, k.param = k.param, annoy.metric = "cosine")
  obj
}

# Function to use Seurat functions to call mitochondrial cluster clones
seuratSNN_cosineDistance <- function(obj, resolution){
  clusters <- FindClusters(object = obj$snn, resolution = resolution)
  return(as.character(clusters[,1]))
}

mut_se_pbmc <- readRDS("../output/filteredpbmcs_mgatk_calls.rds")
mut_se_cd34 <- readRDS("../output/filteredCD34_mgatk_calls.rds")

afin1 <- data.matrix(assays(mut_se_pbmc)[["allele_frequency"]])
afin2 <- data.matrix(assays(mut_se_cd34)[["allele_frequency"]])
afin <- cbind(afin1, afin2)

# Since the PBMCs alone undercalled the mutations, we will examine the mutations called in the PBMCs and CD34s together
obj <- getNN(t(sqrt(afin)))
clusters <- seuratSNN_cosineDistance(obj, resolution = 3.5) 
cluster_name <- str_pad(as.character(clusters), 3, pad = "0")
table(cluster_name)

# Import meta data
sdf <- readRDS("../output/PBMCatac_SignacSeurat_labelTransfer.rds")
sdf$mito_cluster <- cluster_name[1:dim(afin1)[2]]

load("../output/CD34_umap_embedding_granja.rda")
plot_df$mito_cluster <- cluster_name[(dim(afin1)[2]+1): length(cluster_name)]


make_informative_variant_plot <- function(cluster_names, afin, what){
  
  mdf <- data.frame(cluster_names, t(afin)) %>%
    group_by(cluster_names) %>% dplyr::filter(n() >= 5) %>% summarize_all(.funs = mean)
  
  melt_df <- reshape2::melt(mdf, id.vars = "cluster_names")

  # This series of arrangements is for aesthetics
  melt_df2 <- melt_df %>% mutate(g5 = value >= 0.05, g1 = value >= 0.01) %>%  group_by(variable) %>%
    mutate(total_g1 = sum(g1)) %>% 
    arrange(desc((g5)), desc((g1))) %>% dplyr::filter(total_g1 >0)
  
  # More aestetics
  melt_df2$cluster_names <- factor(as.character(melt_df2$cluster_names ), levels = rev(unique( as.character(melt_df2$cluster_names))))
  melt_df2$variable <- factor(as.character(melt_df2$variable ), levels = (unique( as.character(melt_df2$variable))))
  melt_df2$value <- ifelse(melt_df2$value > 0.1, 0.1, melt_df2$value)
  melt_df2$value <- ifelse(melt_df2$value < 0.005, 0, melt_df2$value)
  
  # Visualize heatmap of clustres x variants
  p1 <- ggplot(melt_df2 , aes(x = variable, y = cluster_names, fill = value)) +
    geom_tile() +
    scale_fill_gradientn(colors = jdb_palette("solar_rojos")) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    L_border() +
    theme(legend.position = "none")
  
  cowplot::ggsave(p1, file = paste0("../plots/mutations_clones_",what,"_invivo_grid.pdf"), width = 3.2, height = 1.5)
  print(length(unique(melt_df2$variable)))
  print(length(unique(melt_df2$cluster_names)))
}


make_informative_variant_plot(plot_df$mito_cluster,afin2,"CD34")
make_informative_variant_plot(sdf$mito_cluster,afin1,"PBMC")

saveRDS(plot_df, file = "../output/CD34_clone_DF.rds")
saveRDS(sdf, file = "../output/PBMC_clone_DF.rds")
