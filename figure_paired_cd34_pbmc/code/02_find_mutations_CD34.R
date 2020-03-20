library(SummarizedExperiment)
library(Matrix)
library(dplyr)
library(BuenColors)
library(ComplexHeatmap)
library(circlize)
library(data.table)
library(ggrepel)
source("../../global_functions/variant_calling.R")

"%ni%" <- Negate("%in%")
perc.rank <- function(x) trunc(rank(x))/length(x)

# Import mgatk object
SE <- cbind(readRDS("../../../mtscATACpaper_large_data_files/source/mgatk_output/CD34_G10_v12-mtMask_mgatk.rds"), 
            readRDS("../../../mtscATACpaper_large_data_files/source/mgatk_output/CD34_H8_v12-mtMask_mgatk.rds"))
dim(SE)

# Call variants
mut_se <- call_mutations_mgatk(SE)
misc_df <- data.frame(rowData(mut_se))
vars <- misc_df %>%  filter(n_cells_detected >= 5 & strand_correlation > 0.65 & log10(vmr) > -2 & mean_coverage >= 20) %>% pull(variant)
mut_se_filt <- mut_se[vars,]

# Mark transitions
transition <- c("C>T", "G>A", "A>G", "T>C")
misc_df$transition <- misc_df$nucleotide %in% transition

p1 <- ggplot(misc_df %>%  filter(n_cells_detected >= 5 & mean_coverage >= 20), aes(x = strand_correlation, y = log10(vmr), color = transition)) +
  geom_point(size = 0.4) + scale_color_manual(values = c("black", "dodgerblue")) +
  labs(color = "HQ", x = "Strand concordance", y = "log VMR") +
  pretty_plot(fontsize = 8) + L_border() +
  theme(legend.position = "bottom") + 
  geom_vline(xintercept = 0.65, linetype = 2) +
  geom_hline(yintercept = -2, linetype = 2) + theme(legend.position = "none")
cowplot::ggsave2(p1, file = "../plots/invivoCD34_called_variants.pdf", width = 1.7, height = 1.7)

saveRDS(mut_se_filt, file = "../output/filteredCD34_mgatk_calls.rds")

