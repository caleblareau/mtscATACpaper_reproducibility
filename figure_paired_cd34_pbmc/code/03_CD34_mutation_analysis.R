library(data.table)
library(dplyr)
library(SummarizedExperiment)
library(Matrix)
library(BuenColors)
"%ni%" <- Negate("%in%")

mgatk_se_cd34 <- readRDS("../output/filteredCD34_mgatk_calls.rds")
load("../output/CD34_umap_embedding_granja.rda")

# Find cells in mc9 (HSCs)
hscs <- rownames(plot_df)[plot_df$Clusters == "mc9"]
af_HSCs <- assays(mgatk_se_cd34)[["allele_frequency"]][,colnames(mgatk_se_cd34) %in% hscs]
af_HPCs <- assays(mgatk_se_cd34)[["allele_frequency"]][,colnames(mgatk_se_cd34) %ni% hscs]

# Find mutations enriched in hspcs
df <- data.frame(
  mutation = rownames(af_HSCs),
  af_HSC_mean = rowMeans(af_HSCs), 
  af_HPC_mean = rowMeans(af_HPCs)
) %>% mutate(FC = af_HSC_mean/af_HPC_mean) %>% arrange(desc(FC))

big_plot_df <- data.frame(
  plot_df,
  data.matrix(t(assays(mgatk_se_cd34)[["allele_frequency"]]))
)

ggplot(big_plot_df %>% arrange((X11318T.G)), aes(x = X1, y = X2, color = X11318T.G)) +
  geom_point() + scale_color_gradientn(colors = c("grey", "firebrick")) +
  pretty_plot() + L_border() + theme(legend.position = "bottom")

ggplot(big_plot_df %>% arrange((X1538G.A)), aes(x = X1, y = X2, color = X1538G.A)) +
  geom_point() + scale_color_gradientn(colors = c("grey", "firebrick")) +
  pretty_plot() + L_border() + theme(legend.position = "bottom")

ggplot(big_plot_df %>% arrange((X16131T.C)), aes(x = X1, y = X2, color = X16131T.C)) +
  geom_point() + scale_color_gradientn(colors = c("grey", "firebrick")) +
  pretty_plot() + L_border() + theme(legend.position = "bottom")


p1 <- ggplot(projection_df[dim(projection_df)[1]:1,], aes(x= umap1, y = umap2, color = celltype, label = celltype)) +
  geom_point(size = 0.5) +
  labs(x = "UMAP1", y= "UMAP2", color = "C1 FACS ") +
  pretty_plot() + L_border() + theme(legend.position = "bottom") +
  scale_color_manual(values = c(ejc_color_maps, "none" = "lightgrey", "Monocytes" = "orange2"))
p1

p0 <- ggplot(plot_df %>% filter(Clusters == "mc9"), aes(x= X1, y = X2, color = Clusters)) +
  geom_point(size = 0.5) +
  labs(x = "UMAP1", y= "UMAP2", color = "") +
  pretty_plot() + L_border() + theme(legend.position = "bottom") +
  scale_color_manual(values = jdb_palette('lawhoops'))
p0


