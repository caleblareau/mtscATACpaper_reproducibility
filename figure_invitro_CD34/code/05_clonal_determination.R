library(data.table)
library(Matrix)
library(SummarizedExperiment)
library(Seurat)
library(seriation)
load("trajectory_inferences.17november.rda")

afin <- assays(readRDS("misc_mito/filtered_SE_CD34-500.rds"))[["allele_frequency"]]
seuratSNN <- function(matSVD, resolution = 0.8,k.param = 20){
  set.seed(1)
  rownames(matSVD) <- make.unique(rownames(matSVD))
  obj <- FindNeighbors(matSVD, k.param = k.param)
  clusters <- FindClusters(object = obj$snn, resolution = resolution)
  return(as.character(clusters[,1]))
}

afp <- afin
afp[afp > 0.1] <- 0.1

clusters <- seuratSNN(t(sqrt(afp)), resolution = 2, k.param = 10) 

mdf <- data.frame(clusters, t(afin)) %>%
  group_by(clusters) %>% summarize_all(.funs = mean)

melt_df <- reshape2::melt(mdf, id.vars = "clusters")
melt_df2 <- melt_df %>% mutate(g5 = value >= 0.05, g1 = value >= 0.01) %>%  group_by(variable) %>%
  mutate(total_g1 = sum(g1)) %>% 
  arrange(desc((g5)), desc((g1))) %>% dplyr::filter(total_g1 >0)
melt_df2$clusters <- factor(as.character(melt_df2$clusters ), levels = rev(unique( as.character(melt_df2$clusters))))
melt_df2$variable <- factor(as.character(melt_df2$variable ), levels = (unique( as.character(melt_df2$variable))))
melt_df2$value <- ifelse(melt_df2$value > 0.1, 0.1, melt_df2$value)
melt_df2$value <- ifelse(melt_df2$value < 0.005, 0, melt_df2$value)


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

cowplot::ggsave(p1, file = "plots/i500_grid.pdf", width = 3.2, height = 1.5)

length(unique( as.character(melt_df2$clusters)))
length(unique( as.character(melt_df2$variable)))

saveRDS(data.frame(barcode = colnames(afin), cluster = as.character(clusters)), 
        file = "misc_mito/cluster-id-500.rds")
saveRDS(melt_df2 %>% dplyr::filter(value > 0.005), file = "misc_mito/cluster-mutations-500.rds")
