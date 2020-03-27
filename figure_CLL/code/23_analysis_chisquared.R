library(data.table)
library(dplyr)

# Set up a data frame for plotting
dt <- readRDS("../output/PT1_X2.rds")
dt <- dt %>% arrange(desc(obs_x2)) %>%
  mutate(rank = 1:n(), padj = p.adjust(obs_p), obs_log10p = -1*log10(obs_p))
dt$perm_x2 <- sort(dt$perm_x2, decreasing = TRUE)
dt$perm_log10p <- -1*log10(sort(dt$perm_p, decreasing = FALSE))

sum(dt$padj < 0.01)

p1 <- ggplot(dt, aes(x = rank, y = obs_x2, color = padj < 0.01)) + 
  geom_point_rast(size = 0.1, raster.dpi = 500) + scale_color_manual(values = c("black", "firebrick"))  + 
  geom_point_rast(inherit.aes = FALSE, data = dt, aes(x = rank, y = perm_x2), color = "lightgrey", size = 0.1, raster.dpi = 500) +
  labs(x = "Rank sorted peaks", y = "X2 Statistic") + 
  pretty_plot(fontsize = 8) + L_border() + 
  theme(legend.position = "none")
cowplot::ggsave2(p1, file = "../plots/PT1_ChiSquareSummary.pdf", width = 1.8, height = 1.8)

# Set up a data frame for plotting
dt <- readRDS("../output/PT2_X2.rds")
dt <- dt %>% arrange(desc(obs_x2)) %>%
  mutate(rank = 1:n(), padj = p.adjust(obs_p), obs_log10p = -1*log10(obs_p))
dt$perm_x2 <- sort(dt$perm_x2, decreasing = TRUE)
dt$perm_log10p <- -1*log10(sort(dt$perm_p, decreasing = FALSE))

sum(dt$padj < 0.01)

p1 <- ggplot(dt, aes(x = rank, y = obs_x2, color = padj < 0.01)) + 
  geom_point_rast(size = 0.1, raster.dpi = 500) + scale_color_manual(values = c("black", "firebrick"))  + 
  geom_point_rast(inherit.aes = FALSE, data = dt, aes(x = rank, y = perm_x2), color = "lightgrey", size = 0.1, raster.dpi = 500) +
  labs(x = "Rank sorted peaks", y = "X2 Statistic") + 
  pretty_plot(fontsize = 8) + L_border() + 
  theme(legend.position = "none")
cowplot::ggsave2(p1, file = "../plots/PT2_ChiSquareSummary.pdf", width = 1.8, height = 1.8)