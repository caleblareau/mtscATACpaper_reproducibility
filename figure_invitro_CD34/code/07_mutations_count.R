library(dplyr)
library(SummarizedExperiment)
library(Matrix)
library(BuenColors)

c500 <- readRDS("../output/filtered_mitoSE_CD34-500.rds")
c500_1<- Matrix::colSums(assays(c500)[["allele_frequency"]] >= 0.01)
c500_5<- Matrix::colSums(assays(c500)[["allele_frequency"]] >= 0.05)
c500_10 <- Matrix::colSums(assays(c500)[["allele_frequency"]] >= 0.10)
c500_20 <- Matrix::colSums(assays(c500)[["allele_frequency"]] >= 0.20)

c800 <- readRDS("../output/filtered_mitoSE_CD34-800.rds")
c800_1<- Matrix::colSums(assays(c800)[["allele_frequency"]] >= 0.01)
c800_5<- Matrix::colSums(assays(c800)[["allele_frequency"]] >= 0.05)
c800_10<- Matrix::colSums(assays(c800)[["allele_frequency"]] >= 0.10)
c800_20<- Matrix::colSums(assays(c800)[["allele_frequency"]] >= 0.20)

# function taking a value (proportion of cells) and vector of muts/cell
vv <- function(value, vec){
  sum(vec>=value)/length(vec) * 100
}

df <- rbind(data.frame(
  n1 = sapply(0:10, vv, c500_1),
  n5 = sapply(0:10, vv, c500_5),
  n10 = sapply(0:10, vv, c500_10),
  n20 = sapply(0:10, vv, c500_20),
  what = "500 input (175 mutations)",
  n = 0:10
), 
data.frame(
  n1 = sapply(0:10, vv, c800_1),
  n5 = sapply(0:10, vv, c800_5),
  n10 = sapply(0:10, vv, c800_10),
  n20 = sapply(0:10, vv, c800_20),
  what = "800 input (305 mutations)",
  n = 0:10
))

# Make plot of CDF=like visualization
mdf <- reshape2::melt(df, id.vars = c("what", "n"))
p_out <- ggplot(mdf %>% dplyr::filter(variable != "n20"), aes(x = n, y = value, color = variable)) +
  facet_wrap(~what, nrow = 1) +
  geom_point(size = 0.8) + geom_line() +
  pretty_plot(fontsize = 8) + L_border() +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = c(0:10)) +
  labs(x = "# of mutations per cell", y = "% of cells with mutations") +
  scale_color_manual(values = c("firebrick", "black", "dodgerblue3", "purple3"))

cowplot::ggsave2(p_out, file = "../plots/cells_CDF_like.pdf", width = 3.5, height = 1.8)

# For the callout in the figure
df %>% filter(n == 1)
