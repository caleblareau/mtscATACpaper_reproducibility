library(data.table)
library(dplyr)

pt1_vars <- readRDS("output/PT1_specificVariants_forSupplement.rds")
pt2_vars <- readRDS("output/PT2_specificVariants_forSupplement.rds")

plot_df <- data.frame(
  het = c(pt1_vars$X5140G.A, pt1_vars$X14858G.A,
          pt1_vars$X1872T.C, pt1_vars$X1260A.G,
          pt2_vars$X12980G.A, pt2_vars$X4853G.A) * 100,
  mut = c(rep("5140G>A", dim(pt1_vars)[1]), rep("14858G>A", dim(pt1_vars)[1]), 
          rep("1872T>C", dim(pt1_vars)[1]), rep("1260A>G", dim(pt1_vars)[1]), 
          rep("12980G>A", dim(pt2_vars)[1]), rep("4853G>A", dim(pt2_vars)[1]))
)

plot_df$mut <- factor(as.character(plot_df$mut), levels = c("5140G>A", "14858G>A", "4853G>A", "1872T>C", "1260A>G", "12980G>A"))
p1 <- ggplot(plot_df, aes(x = het)) +
  geom_histogram(binwidth = 10, fill = "black") + facet_wrap(~mut ) +
  pretty_plot(fontsize = 7) 
cowplot::ggsave(p1, file = "output/mut_histograms.pdf",width = 3.5, height = 1.8)
