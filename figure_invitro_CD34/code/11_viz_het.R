library(data.table)
library(dplyr)
library(ggrastr)
load("../output/trajectory_inferences.18march2020.rda")

mut_se <- readRDS(paste0("../output/filtered_mitoSE_CD34-800.rds"))


muts <- c(
  "12451A>G", "3022G>A", "791G>A",  "11682G>C", "3047G>A"
)

df <- data.frame(
  colData(mut_se), data.matrix(t(assays(mut_se)[["allele_frequency"]][muts,])))
mdf <- merge(df, plot_df, by = "row.names")

ggplot(mdf %>% arrange(X3047G.A), aes(x = X1, y = X2, color = X3047G.A > 0.1)) +
  geom_point() + facet_wrap(~ library) + scale_color_manual(values = c("lightgrey", "firebrick")) +
  pretty_plot() +  theme(legend.position = "bottom") + labs(x = "UMAP1", y = "UMAP2")



ggplot(mdf %>% arrange(X10510T.C), aes(x = X1, y = X2, color = X10510T.C > 0.01)) +
  geom_point() + facet_wrap(~ library) + scale_color_manual(values = c("lightgrey", "firebrick")) +
  pretty_plot() +  theme(legend.position = "bottom") + labs(x = "UMAP1", y = "UMAP2")

ggplot(mdf, aes(x = X12962G.A, y = X10510T.C)) +
  geom_point()


ggplot(mdf %>% arrange(X12962G.C), aes(x = X1, y = X2, color = X12962G.C > 0.01)) +
  geom_point() + facet_wrap(~ library) + scale_color_manual(values = c("lightgrey", "firebrick")) +
  pretty_plot() +  theme(legend.position = "bottom") + labs(x = "UMAP1", y = "UMAP2")



ggplot(mdf %>% arrange(X791G.A), aes(x = X1, y = X2, color = X791G.A > 0.05)) +
  geom_point() + facet_wrap(~ library) + scale_color_manual(values = c("lightgrey", "firebrick")) +
  pretty_plot() +  theme(legend.position = "bottom") + labs(x = "UMAP1", y = "UMAP2")

ggplot(mdf %>% arrange(X3022G.A), aes(x = X1, y = X2, color = X3022G.A > 0.05)) +
  geom_point() + facet_wrap(~ library) + scale_color_manual(values = c("lightgrey", "firebrick")) +
  pretty_plot() +  theme(legend.position = "bottom") + labs(x = "UMAP1", y = "UMAP2")

ggplot(mdf %>% arrange(X12451A.G), aes(x = X1, y = X2, color = X12451A.G > 0.05)) +
  geom_point() + facet_wrap(~ library) + scale_color_manual(values = c("lightgrey", "firebrick")) +
  pretty_plot() +  theme(legend.position = "bottom") + labs(x = "UMAP1", y = "UMAP2")


ggplot(mdf, aes(x = X3022G.A, y = X791G.A, color = X12451A.G > 0.01)) +
  geom_point()


