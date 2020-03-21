library(dplyr)
library(data.table)
library(BuenColors)

import_complexity <- function(file, description){
  dt <- data.frame(fread(file))
  dt$library <- description
  dt %>% arrange(ESTIMATED_LIBRARY_SIZE)
}

top1k_df <- rbind(
  import_complexity("../data/mgatk_complexity/Mix_1_CR-mtMask_mgatk_mitoComplexity.tsv.gz", "aRegular"),
  import_complexity("../data/mgatk_complexity/Mix_2_NP40_CR-mtMask_mgatk_mitoComplexity.tsv.gz", "bNoFA"),
  import_complexity("../data/mgatk_complexity/Mix_Fix_1h_CR-mtMask_mgatk_mitoComplexity.tsv.gz", "FA"),
  import_complexity("../data/mgatk_complexity/Mix_Fix_6h_CR-mtMask_mgatk_mitoComplexity.tsv.gz", "FA")
) %>% group_by(library) %>% top_n(1000, wt = ESTIMATED_LIBRARY_SIZE)


top1k_df %>% group_by(library) %>% summarize(median(ESTIMATED_LIBRARY_SIZE))

pmtComp <- ggplot(top1k_df, aes(x = library, y = log10(ESTIMATED_LIBRARY_SIZE))) +
  geom_boxplot(outlier.shape = NA, width = 0.8) +
  coord_cartesian(ylim = c(3.0, 5.5))  + scale_y_continuous(expand = c(0,0)) + 
  pretty_plot(fontsize = 6) + theme(legend.position = "bottom") + labs(x = "", color = "", y = "log10 mtDNA complexity") +
  scale_color_manual(values = "black") + L_border()+ theme(axis.text.x = element_text(angle = 0, hjust = 1))

cowplot::ggsave2(pmtComp, file = "../plots/mtDNA_complexity.pdf", width = 1.8, height = 1.8)
