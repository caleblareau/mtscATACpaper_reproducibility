library(data.table)
library(Matrix)
library(SummarizedExperiment)
library(ggbeeswarm)
library(stringr)
load("../output/trajectory_inferences.18march2020.rda")

afin <- assays(readRDS("../output/filtered_mitoSE_CD34-800.rds"))[["allele_frequency"]]

# Subset to variants for visualization
plot_vars <- c("3243A>T", "12316G>A")
plot_vars_ugly <- paste0("X", gsub(">", ".", plot_vars))
plotraw_df <- reshape2::melt(data.matrix(afin[plot_vars,])) %>% data.frame() 

# Visualize cells with at least 2% heteroplasmy
pdf <- reshape2::melt(data.matrix(afin[plot_vars,])) %>% data.frame() %>% dplyr::filter(value >= 0.02) %>%
  mutate(blah = str_split_fixed(Var2, "-", 3)[,1]) %>% mutate(timepoint = factor(as.numeric(factor(blah))))

# Make plot
p1 <- ggplot(pdf, aes(x = timepoint, y = value*100, color = timepoint)) +
  facet_wrap(~Var1, scales = "free_y", nrow = 1) +
  geom_quasirandom(size = 0.1) +
  scale_color_manual(values = c("dodgerblue", "#008080", "purple3")) +
  pretty_plot(fontsize = 8) + L_border() +
  labs(x = "Timepoint", y = "Heteroplasmy") +
  theme(legend.position = "none")
cowplot::ggsave2(p1, file = "../plots/CD34_800_quasirandom_new.pdf", width = 2.2, height = 1.7)

