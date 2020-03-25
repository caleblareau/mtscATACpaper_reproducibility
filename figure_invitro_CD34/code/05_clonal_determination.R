library(data.table)
library(Matrix)
library(SummarizedExperiment)
library(Seurat)
library(seriation)
load("../output/trajectory_inferences.18march2020.rda")

# Function to use Seurat functions to call mitochondrial cluster clones
seuratSNN_cosineDistance <- function(mat_af, resolution = 0.8,k.param = 20){
  set.seed(1)
  rownames(mat_af) <- make.unique(rownames(mat_af))
  obj <- FindNeighbors(mat_af, k.param = k.param, annoy.metric = "cosine")
  clusters <- FindClusters(object = obj$snn, resolution = resolution)
  return(as.character(clusters[,1]))
}

# Master function to determine the clone assignments for each 
process_colonies <- function(cell_input_n){
  
  # Import allele frequencies
  afin <- assays(readRDS(paste0("../output/filtered_mitoSE_CD34-",cell_input_n,".rds")))[["allele_frequency"]]
  
  # Call clone/clusters
  clusters <- seuratSNN_cosineDistance(t(sqrt(afin)), resolution = 1.5, k.param = 10) 
  
  # Kind of nasty data frame but the idea here is to build and order mutations
  # In clones in order
  mdf <- data.frame(clusters, t(afin)) %>%
    group_by(clusters) %>% summarize_all(.funs = mean)
  
  melt_df <- reshape2::melt(mdf, id.vars = "clusters")
  
  # This series of arrangements is for aesthetics
  melt_df2 <- melt_df %>% mutate(g5 = value >= 0.05, g1 = value >= 0.01) %>%  group_by(variable) %>%
    mutate(total_g1 = sum(g1)) %>% 
    arrange(desc((g5)), desc((g1))) %>% dplyr::filter(total_g1 >0)
  
  # More aestetics
  melt_df2$clusters <- factor(as.character(melt_df2$clusters ), levels = rev(unique( as.character(melt_df2$clusters))))
  melt_df2$variable <- factor(as.character(melt_df2$variable ), levels = (unique( as.character(melt_df2$variable))))
  melt_df2$value <- ifelse(melt_df2$value > 0.1, 0.1, melt_df2$value)
  melt_df2$value <- ifelse(melt_df2$value < 0.005, 0, melt_df2$value)
  
  # Visualize heatmap of clustres x variants
  p1 <- ggplot(melt_df2, aes(x = variable, y = clusters, fill = value)) +
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
  
  cowplot::ggsave(p1, file = paste0("../plots/mutations_clones_",cell_input_n,"_grid.pdf"), width = 3.2, height = 1.5)
  
  # Export assignments for downstream use
  saveRDS(data.frame(barcode = colnames(afin), cluster = as.character(clusters)), 
          file = paste0("../output/cluster-id-",cell_input_n,".rds"))
  saveRDS(melt_df2 %>% dplyr::filter(value > 0.005), file = paste0("../output/cluster-mutations-",cell_input_n,".rds"))
  
  # Report the number of clones called and the number of informative variants used in the assignment
  c(length(unique( as.character(melt_df2$clusters))), length(unique( as.character(melt_df2$variable))))
}
process_colonies("500")
process_colonies("800")
