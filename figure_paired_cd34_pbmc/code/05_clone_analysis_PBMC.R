library(data.table)
library(dplyr)
library(SummarizedExperiment)
library(BuenColors)

# Function to use Seurat functions to call mitochondrial cluster clones
seuratSNN_cosineDistance <- function(mat_af, resolution = 0.8,k.param = 20){
  set.seed(1)
  rownames(mat_af) <- make.unique(rownames(mat_af))
  obj <- FindNeighbors(mat_af, k.param = k.param, annoy.metric = "cosine")
  clusters <- FindClusters(object = obj$snn, resolution = resolution)
  return(as.character(clusters[,1]))
}

mut_se <- readRDS("../output/filteredpbmcs_mgatk_calls.rds")
sdf <- readRDS("../output/PBMCatac_SignacSeurat_labelTransfer.rds")

afin <- data.matrix(assays(mut_se)[["allele_frequency"]])
clusters <- seuratSNN_cosineDistance(t(sqrt(afin)), resolution = 2, k.param = 10) 
table(clusters)

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

cowplot::ggsave(p1, file = paste0("../plots/mutations_clones_PBMCs_invivo_grid.pdf"), width = 3.2, height = 1.5)

for_p <- data.frame(clusters, sdf) %>%
  dplyr::filter(lineage_assignment_atac != "other") %>% 
  group_by(clusters, lineage_assignment_atac) %>% summarize(count = n()) %>%
  ungroup() %>% group_by(clusters) %>% mutate(total = sum(count)) %>%
  mutate(proportion = count/total) %>%
  dplyr::filter(total > 10) %>% arrange(desc(proportion))

for_p$ca <- factor(as.character(for_p$clusters), levels = unique(as.character(for_p$clusters)))

ggplot(for_p, aes(fill=lineage_assignment_atac, y=total, x=ca)) + 
  geom_bar(position="stack", stat="identity")

ggplot(for_p, aes(fill=lineage_assignment_atac, y=proportion, x=ca)) + 
  geom_bar(position="stack", stat="identity")

mito_prop_df <- data.frame(clusters, sdf, colData(mut_se))

ggplot(mito_prop_df, aes(x = predicted.id, y = depth/nCount_ATAC *70)) + 
  geom_boxplot()


