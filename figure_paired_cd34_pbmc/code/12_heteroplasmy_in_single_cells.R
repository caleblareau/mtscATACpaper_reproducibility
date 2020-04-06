library(Matrix)
library(SummarizedExperiment)
library(BuenColors)
library(dplyr)

# Import heteroplasmy data

mgatk_se_cd34 <- readRDS("../output/filteredCD34_mgatk_calls.rds")
mgatk_se_pbmc <- readRDS("../output/filteredpbmcs_mgatk_calls.rds")

# Find the cells/variants with minimum 20x coverage
rm_cd34 <- rowMaxs(data.matrix((assays(mgatk_se_cd34)[["coverage"]] > 20) * (assays(mgatk_se_cd34)[["allele_frequency"]])))
rm_pbmc <- rowMaxs(data.matrix((assays(mgatk_se_pbmc)[["coverage"]] > 20) * (assays(mgatk_se_pbmc)[["allele_frequency"]])))


p1 <- ggplot(data.frame(af=rm_cd34), aes(x = af*100)) +
  geom_histogram(bins = 20, fill = "lightgrey", color = "black") +
  pretty_plot(fontsize = 7) + L_border() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  ggtitle("CD34")

p2 <- ggplot(data.frame(af=rm_pbmc), aes(x = af*100)) +
  geom_histogram(bins = 20, fill = "lightgrey", color = "black") +
  pretty_plot(fontsize = 7) + L_border() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  ggtitle("PBMC")

cowplot::ggsave2(cowplot::plot_grid(p1, p2, nrow = 2), width = 1.6, height = 2, filename = "../plots/singlecell_heteroplasmy_histograms.pdf")

sum(rm_pbmc > 0.95)
sum(rm_cd34 > 0.95)
