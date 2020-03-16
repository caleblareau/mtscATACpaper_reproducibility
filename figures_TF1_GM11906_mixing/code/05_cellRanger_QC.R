library(dplyr)
library(data.table)
library(BuenColors)
library(ggrastr)
library(SummarizedExperiment)
library(Matrix)

source("../../global_functions/estimateLibrarySize.R")

# Get top n cell qc metrics per run
importQC <- function(raw_name, name, n = 1000){
  
  dt <- fread(paste0("../data/cellrangeratac_singlecell/",raw_name,"_singlecell.csv.gz")) %>% 
    dplyr::filter(cell_id != "None") %>% 
    mutate(Experiment = name) %>%
    mutate(DNaseProp = DNase_sensitive_region_fragments/passed_filters) %>%
    mutate(MitoProp = mitochondrial/total) %>% 
    mutate(duplicateProp = duplicate / total) %>%
    mutate(TSSProp = TSS_fragments/passed_filters) %>% data.frame()
  dt$LibraryComplexity <- sapply(1:dim(dt)[1], function(i){
    estimateLibrarySize(dt[i,"total"],dt[i,"passed_filters"])})
  
  top_1k <- dt %>% top_n(n, LibraryComplexity)
  top_1k$barcode <- as.character(top_1k$barcode)
  return(top_1k)

}

screen_df <- rbindlist(list(
  importQC("Mix_1", "aRegular"),
  importQC("Mix_2_NP40", "Condition1"),
  importQC("Mix_3_Fix_1h", "Condition2"),
  importQC("Mix_3_Fix_12h", "Condition3"),
  importQC("Mix_4_Fix_1h", "Condition4"),
  importQC("Mix_4_Fix_12h", "Condition5")
))

# Summarize values for the initial screen of conditions

screen_df %>% group_by(Experiment) %>%
  summarize(median(MitoProp)*100, median(DNaseProp)*100, median(TSSProp)*100, median(LibraryComplexity))

pA <- ggplot(screen_df, aes(x = Experiment, y = MitoProp*100)) +
  geom_boxplot(outlier.shape = NA, fill = NA, width = 0.2) + geom_violin(fill = NA) +
  coord_cartesian(ylim = c(0, 50)) + scale_y_continuous(expand = c(0,0)) + 
  pretty_plot(fontsize = 7) + theme(legend.position = "bottom") + labs(x = "", color = "", y = "% mtDNA fragments") +
  scale_color_manual(values = "black") + L_border() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

pB <- ggplot(screen_df, aes(x = Experiment, y = DNaseProp*100)) +
  geom_boxplot(outlier.shape = NA, fill = NA, width = 0.2) + geom_violin(fill = NA) +
  coord_cartesian(ylim = c(35, 90)) + scale_y_continuous(expand = c(0,0)) + 
  pretty_plot(fontsize = 7) + theme(legend.position = "bottom") + labs(x = "", color = "", y = "% reads in DNase") +
  scale_color_manual(values = "black") + L_border() + theme(axis.text.x = element_text(angle = 45, hjust = 1))


cowplot::ggsave2(cowplot::plot_grid(pA, pB, nrow = 1), file = "../plots/Panel1BC.pdf", width = 3.7, height = 2)

# Summarize additional plots for the supplement

pC <- ggplot(screen_df, aes(x = Experiment, y = TSSProp*100)) +
  geom_boxplot(outlier.shape = NA, fill = NA, width = 0.2) + geom_violin(fill = NA) +
  coord_cartesian(ylim = c(0, 65)) + scale_y_continuous(expand = c(0,0)) + 
  pretty_plot(fontsize = 7) + theme(legend.position = "bottom") + labs(x = "", color = "", y = "% reads in TSS") +
  scale_color_manual(values = "black") + L_border() + theme(axis.text.x = element_text(angle = 0, hjust = 1))

cowplot::ggsave2(pC, file = "../plots/tss_rates.pdf", width = 1.8, height = 1.8)

pD <-  ggplot(screen_df, aes(x = Experiment, y = log10(LibraryComplexity))) +
  geom_boxplot(outlier.shape = NA, fill = NA, width = 0.2) + geom_violin(fill = NA) +
  coord_cartesian(ylim = c(3.5, 5.5))  + scale_y_continuous(expand = c(0,0)) + 
  pretty_plot(fontsize = 6) + theme(legend.position = "bottom") + labs(x = "", color = "", y = "log10 library size") +
  scale_color_manual(values = "black") + L_border()+ theme(axis.text.x = element_text(angle = 0, hjust = 1))

cowplot::ggsave2(pD, file = "../plots/library_complexity.pdf", width = 1.8, height = 1.8)

# Characterize the effects of mtDNA masking

mtMappingRatesDF <- rbindlist(list(
  importQC("Mix_1", "aRegular_NoMask"), 
  importQC("Mix_1_CR-mtMask", "aRegular_Mask"),
  importQC("Mix_2_NP40", "bNoFA_NoMask"),
  importQC("Mix_2_NP40_CR-mtMask", "bNoFA_Mask"),
  importQC("Mix_Fix_1h", "cwFA_NoMask", n = 500),
  importQC("Mix_Fix_1h_CR-mtMask", "cwFA_Mask", n = 500),
  importQC("Mix_Fix_6h", "cwFA_NoMask", n = 500),
  importQC("Mix_Fix_6h_CR-mtMask", "cwFA_Mask", n = 500)
))

mtMappingRatesDF %>% group_by(Experiment) %>%
  summarize(mDNase = median(DNaseProp), mMitoProp = median(MitoProp), mLibraryComplexity = median(LibraryComplexity))

two <- stringr::str_split_fixed(as.character(mtMappingRatesDF$Experiment), "_", 2)
mtMappingRatesDF$library <- two[,1]
mtMappingRatesDF$mask <- factor(two[,2], levels = c("NoMask", "Mask"))

pLC <- ggplot(mtMappingRatesDF, aes(x = library, y = log10(LibraryComplexity), color = mask)) +
  geom_boxplot(outlier.shape = NA, width = 0.8) +
  scale_y_continuous(limits = c(3.5, 5.5), expand = c(0,0)) +
  pretty_plot(fontsize = 8) + L_border() +
  labs(x = "", y = "log10 chromatin complexity") +
  scale_color_manual(values = c("dodgerblue3", "firebrick"))

pMR <- ggplot(mtMappingRatesDF, aes(x = library, y = MitoProp*100, color = mask)) +
  geom_boxplot(outlier.shape = NA, width = 0.8) +
  pretty_plot(fontsize = 8) + L_border() +
  scale_y_continuous( expand = c(0,0)) +
  labs(x = "", y = "% reads aligned to mtDNA") +
  scale_color_manual(values = c("dodgerblue3", "firebrick"))

two_plot <- cowplot::ggsave2(cowplot::plot_grid(pMR, pLC, nrow =1), filename = "../plots/compareMasking.pdf", height = 2, width = 8)

p1 <- ggplot(mtMappingRatesDF %>% filter(mask == "NoMask"), aes(x = library, y = log10(LibraryComplexity))) +
  geom_boxplot(outlier.shape = NA, width = 0.8) +
  scale_y_continuous(limits = c(3.5, 5.5), expand = c(0,0)) +
  pretty_plot(fontsize = 8) + L_border() +
  labs(x = "", y = "log10 chromatin complexity")

cowplot::ggsave2(p1, file = "../plots/FA_library_complexity.pdf", width = 1.8, height = 1.8)



