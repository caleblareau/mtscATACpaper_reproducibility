library(dplyr)
library(data.table)
library(BuenColors)
library(RcppRoll)
library(Matrix)
library(SummarizedExperiment)

# Process the before mtscATAC-seq
SE <- readRDS("../../../mtscATACpaper_large_data_files/source/mgatk_output/Mix_1_mgatk.rds")
top2k <- tail(head(sort(colData(SE)$depth, decreasing = TRUE),2000),1)
cov_original <- rowMeans((assays(SE)[["coverage"]][,top2k < colData(SE)$depth]))
mean(cov_original)

# Import the after
a <- readRDS("../../../mtscATACpaper_large_data_files/source/mgatk_output/Mix_Fix_1h_CR-mtMask_mgatk.rds")
b <- readRDS("../../../mtscATACpaper_large_data_files/source/mgatk_output/Mix_Fix_6h_CR-mtMask_mgatk.rds")
mcols(a)$refAllele <- toupper(as.character(mcols(a)$refAllele))
mcols(b)$refAllele <- toupper(as.character(mcols(b)$refAllele))
SE <- cbind(a,b)

# Process the after
top2k <- tail(head(sort(colData(SE)$depth, decreasing = TRUE),2000),1)
cov_new <- rowMeans((assays(SE)[["coverage"]][,top2k < colData(SE)$depth]))
mean(cov_new)

# Theme to remove any of the other plot riff raff
xxtheme <-   theme(
  axis.line = element_blank(),
  axis.ticks.y = element_blank(),   
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.margin = unit(c(0, 0, 0, 0), "cm"),
  plot.margin = unit(c(-0.35, -0.35, -0.35, -0.35), "cm")) 

# Gental smooth of just 5 bp for plot aesthetics
smooth <- 5
data.frame(
  pos = roll_mean(1:length(cov_new), smooth),
  new = roll_mean(cov_new, smooth),  
  original = roll_mean(cov_original, smooth)
) -> df

mdf <- reshape2::melt(df, id.var = "pos")

# Visualize the rolled means
P1 <- ggplot(mdf, aes(x = pos, y = value, color = variable)) + 
  geom_line() +  expand_limits(y = c(-5, 4)) +
  pretty_plot(fontsize = 8)  + scale_color_manual(values = c("firebrick", "dodgerblue4")) +
  coord_polar(direction = 1) + labs(x = "", y = "log2 Coverage") + scale_y_log10() +
  theme(legend.position = "none") + xxtheme 

cowplot::ggsave2(P1, file = "../plots/rollMean_coverage.pdf", width = 4, height = 4)


# R2R comment-- see what residuate factors affect dips in coverage

# First, compute GC per bin
df2 <- data.frame(
  pos = 1:length(cov_new),
  rowData(SE), 
  cov_new)
df2$isGC <- df2$refAllele %in% c("G", "C")

# Visualize effects of GC
smooth_wide <- 50
roll_mean_df_gc <- data.frame(
  pos = roll_mean(1:length(cov_new), smooth_wide, by = smooth_wide/2),
  coverage = roll_mean(cov_new, smooth_wide, by = smooth_wide/2),
  GC = roll_mean(df2$isGC, smooth_wide, by = smooth_wide/2)
)

cor(roll_mean_df_gc$GC, roll_mean_df_gc$coverage)

p1 <- ggplot(roll_mean_df_gc, aes(x = GC *100, y = coverage)) +
  geom_point(alpha = 0.5, size = 0.5) +
  pretty_plot(fontsize = 8) + L_border() + labs(x = "GC Content (%)", y = "mtDNA coverage") +
  geom_smooth(method='lm', formula= y~x, color = "firebrick", linetype = 2, se = FALSE)
cowplot::ggsave2(p1, file = "../plots/residual_variation_coverage.pdf", width = 1.6, height = 1.6)
