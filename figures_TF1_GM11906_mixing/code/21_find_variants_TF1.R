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
SE <- cbind(readRDS("../../../mtscATACpaper_large_data_files/source/mgatk_output/Mix_Fix_1h_CR-mtMask_mgatk.rds"), 
            readRDS("../../../mtscATACpaper_large_data_files/source/mgatk_output/Mix_Fix_6h_CR-mtMask_mgatk.rds"))

# Find TF1 cells
rbind(read.table("../output/data1_meta.tsv", header = TRUE), 
      read.table("../output/data6_meta.tsv", header = TRUE)) %>% filter(assign == "TF1") %>%
  filter(mean_cov > 50) -> TF_cells_df

SE <- SE[,as.character(TF_cells_df$cell_id)]

# Call variants
mut_se <- call_mutations_mgatk(SE)
misc_df <- data.frame(rowData(mut_se))
filter_df <- misc_df %>%  filter(n_cells_detected >= 5 & strand_correlation > 0.65 & log10(vmr) > -2)
dim(filter_df)

p1 <- ggplot(misc_df %>%  filter(n_cells_detected >= 5 ), aes(x = strand_correlation, y = log10(vmr), color = log10(vmr) > -2 & strand_correlation > 0.65)) +
  geom_point(size = 0.4) + scale_color_manual(values = c("black", "firebrick")) +
  labs(color = "HQ", x = "Strand concordance", y = "log VMR") +
  pretty_plot(fontsize = 8) + L_border() +
  theme(legend.position = "bottom") + 
  geom_vline(xintercept = 0.65, linetype = 2) +
  geom_hline(yintercept = -2, linetype = 2) + theme(legend.position = "none")
cowplot::ggsave2(p1, file = "../plots/TF1_final_var_call_ndetect.pdf", width = 1.7, height = 1.7)

write.table(data.frame(filter_df %>% arrange(desc(strand_correlation)) %>% pull(variant)), 
            file = "../output/TF1_VMR_processed_variants_ndetect.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


