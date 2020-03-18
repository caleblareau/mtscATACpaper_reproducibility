library(data.table)
library(Matrix)
library(SummarizedExperiment)
library(ggbeeswarm)
library(stringr)
load("trajectory_inferences.17november.rda")

afin <- assays(readRDS("misc_mito/filtered_SE_CD34-800.rds"))[["allele_frequency"]]

#plot_vars <- c("3243A>T", "13443T>C", "3712G>A", "14569G>A")
plot_vars <- c("3243A>T", "12316G>A")

plot_vars_ugly <- paste0("X", gsub(">", ".", plot_vars))

plotraw_df <- reshape2::melt(data.matrix(afin[plot_vars,])) %>% data.frame() 

pdf <- reshape2::melt(data.matrix(afin[plot_vars,])) %>% data.frame() %>% dplyr::filter(value >= 0.02) %>%
  mutate(blah = str_split_fixed(Var2, "-", 3)[,1]) %>% mutate(timepoint = factor(as.numeric(factor(blah))))

p1 <- ggplot(pdf, aes(x = timepoint, y = value*100, color = timepoint)) +
  facet_wrap(~Var1, scales = "free_y", nrow = 1) +
  geom_quasirandom(size = 0.1) +
  scale_color_manual(values = c("dodgerblue", "#008080", "purple3")) +
  pretty_plot(fontsize = 8) + L_border() +
  labs(x = "Timepoint", y = "Heteroplasmy") +
  theme(legend.position = "none")
cowplot::ggsave(p1, file = "plots/CD34_800_quasirandom_new.pdf", width = 3.5, height = 1.7)

