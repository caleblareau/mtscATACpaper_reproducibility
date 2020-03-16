library(dplyr)
library(data.table)
library(BuenColors)
library(RcppRoll)
library(Matrix)
library(SummarizedExperiment)

SE <- readRDS("../data/big_gi/Mix_1_mgatk.rds")
top2k <- tail(head(sort(colData(SE)$depth, decreasing = TRUE),2000),1)
cov_original <- rowMeans((assays(SE)[["coverage"]][,top2k < colData(SE)$depth]))
mean(cov_original)

a <- readRDS("../data/big_gi/Mix_Fix_1h_CR-mtMask_mgatk.rds")
b <- readRDS("../data/big_gi/Mix_Fix_6h_CR-mtMask_mgatk.rds")
mcols(a)$refAllele <- toupper(as.character(mcols(a)$refAllele))
mcols(b)$refAllele <- toupper(as.character(mcols(b)$refAllele))
SE <- cbind(a,b)

top2k <- tail(head(sort(colData(SE)$depth, decreasing = TRUE),2000),1)
cov_new <- rowMeans((assays(SE)[["coverage"]][,top2k < colData(SE)$depth]))
mean(cov_new)


xxtheme <-   theme(
  axis.line = element_blank(),
  axis.ticks.y = element_blank(),        ## <- this line
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.margin = unit(c(0, 0, 0, 0), "cm"),
  plot.margin = unit(c(-0.35, -0.35, -0.35, -0.35), "cm")) 

smooth <- 5
data.frame(
  pos = roll_mean(1:length(cov_new), smooth),
  new = roll_mean(cov_new, smooth),  
  original = roll_mean(cov_original, smooth)
) -> df

mdf <- reshape2::melt(df, id.var = "pos")

ggplot(mdf, aes(x = pos, y = value, color = variable)) + 
  geom_line() +  expand_limits(y = c(-5, 4)) +
  pretty_plot(fontsize = 8)  + scale_color_manual(values = c("firebrick", "dodgerblue4")) +
  coord_polar(direction = 1) + labs(x = "", y = "log2 Coverage") + scale_y_log10() +
   theme(legend.position = "none") + xxtheme -> P1

cowplot::ggsave(P1, file = "../plots/rollMean_coverage.pdf", width = 4, height = 4)

