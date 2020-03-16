library(dplyr)
library(data.table)
library(BuenColors)
library(ggrastr)
library(SummarizedExperiment)
library(Matrix)

source("../../global_functions/estimateLibrarySize.R")

merrf <- "8344A>G"
tf1_vars <- c("9824T>C", "6680T>C", "8701A>G", "15301G>A", "12811T>C", "14783T>C", "4164A>G", "7853G>A", 
              "7684T>C", "10345T>C", "12705C>T", "15043G>A", "5351A>G", "10873T>C", "150C>T", "6455C>T",
              "9540T>C", "12405C>T", "199T>C", "16129G>A", "10400C>T", "8552T>C", "4071C>T", "4048G>A",
              "16298T>C", "16297T>C", "16223C>T", "5460G>A", "16189T>C")
gm11906_vars <- c("6734G>A", "3010G>A", "6554C>T", "14798T>C", "11251A>G", "15452C>A", "12612A>G", "16069C>T",
                  "16092T>C", "16126T>C", "12127G>A", "462C>T", "4216T>C", "16261C>T", "228G>A", "295C>T")

assign_pull <- function(SE){
  cov <- assays(SE)[["coverage"]]+ 0.001
  ref_allele <- as.character(rowRanges(SE)$refAllele)
  
  getMutMatrix <- function(letter){
    mat <- (assays(SE)[[paste0(letter, "_counts_fw")]] + assays(SE)[[paste0(letter, "_counts_rev")]]) / cov
    rownames(mat) <- paste0(as.character(1:dim(mat)[1]), toupper(ref_allele), ">", letter)
    return(mat[toupper(ref_allele) != letter,])
  }
  
  mmat <- rbind(getMutMatrix("A"), getMutMatrix("C"), getMutMatrix("G"), getMutMatrix("T"))
  df <- data.frame(
    cell_id = colnames(mmat),
    tf1_mean = round(Matrix::colMeans(mmat[tf1_vars,]), 3),
    gm11906_mean = round(Matrix::colMeans(mmat[gm11906_vars,]), 3),
    merrf_af  = round(mmat[merrf,], 3),
    merrf_cov = round(cov[8344,], 0),
    mean_cov = round(Matrix::colMeans(cov), 1)
  )
  df$assign <- ifelse(
    df$gm11906_mean > 0.90, "GM11906", 
    ifelse(df$tf1_mean > 0.90, "TF1", 
           ifelse(df$mean_cov < 3, "Low_coverage", "Collision")))
  df
  
}

# Import all barcodes passing knee + assign a cell type
importQC <- function(directory, name, n = 1000){
  
  #SE <- readRDS(paste0(directory, "_CR-mtMask_mgatk/final/",directory,"_CR-mtMask_mgatk.rds"))
  #mito_df <- assign_pull(SE)
  
  barcodes <- fread(paste0(directory, "/outs/filtered_peak_bc_matrix/barcodes.tsv"))[[1]]
  dt <- fread(paste0(directory, "/outs/singlecell.csv")) %>% filter(barcode %in% barcodes) %>%
    mutate(Experiment = name) %>%
    mutate(DNaseProp = DNase_sensitive_region_fragments/passed_filters) %>%
    mutate(MitoProp = mitochondrial/total) %>% 
    mutate(duplicateProp = duplicate / total) %>%
    mutate(TSSProp = TSS_fragments/passed_filters) %>% data.frame()
  dt$LibraryComplexity <- sapply(1:dim(dt)[1], function(i){
    estimateLibrarySize(dt[i,"total"],dt[i,"passed_filters"])})
  
  top_1k <- dt %>% top_n(n, LibraryComplexity)
  top_1k$barcode <- as.character(top_1k$barcode)
  #mito_df$cell_id <- as.character(mito_df$cell_id)
  return(top_1k)
  #dd <- merge(top_1k, mito_df, by.x = "barcode", by.y = "cell_id")
  #sreturn(dd)
}

importQCquick <- function(directory, name, n = 1000){
  
  dt <- fread(paste0(directory, "/outs/singlecell.csv")) %>% f
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

screen_df <- compareTwoDF <- rbindlist(list(
  importQCquick("Mix_1", "aRegular"),
  importQCquick("Mix_2_NP40", "Condition1"),
  importQCquick("Mix_3_Fix_1h", "Condition2"),
  importQCquick("Mix_3_Fix_12h", "Condition3"),
  importQCquick("Mix_4_Fix_1h", "Condition4"),
  importQCquick("Mix_4_Fix_12h", "Condition5")
))

screen_df %>% group_by(Experiment) %>%
  summarize(median(MitoProp)*100, median(DNaseProp)*100, median(TSSProp)*100)


ggplot(screen_df, aes(x = Experiment, y = log10(LibraryComplexity))) +
  geom_violin() + geom_boxplot(outlier.shape = NA, fill = NA, width = 0.2) +
  coord_cartesian(ylim = c(3, 6))  + scale_y_continuous(expand = c(0,0)) + 
  pretty_plot(fontsize = 6) + theme(legend.position = "bottom") + labs(x = "", color = "", y = "log10 library size") +
  scale_color_manual(values = "black") + L_border()+ theme(axis.text.x = element_text(angle = 45, hjust = 1))


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



mtMappingRatesDF <- rbindlist(list(
  importQCquick("Mix_1", "aRegular_NoMask"), 
  importQCquick("Mix_1_CR-mtMask", "aRegular_Mask"),
  importQCquick("Mix_2_NP40", "bNoFA_NoMask"),
  importQCquick("Mix_2_NP40_CR-mtMask", "bNoFA_Mask"),
  importQCquick("../merged_reseq/Mix_Fix_1h_CR_merged", "cwFA_NoMask", n = 500),
  importQCquick("../merged_reseq/Mix_Fix_1h_CR_merged-mtMask", "cwFA_Mask", n = 500),
  importQCquick("../17November2018/Mix_Fix_6h_CR/", "cwFA_NoMask", n = 500),
  importQCquick("../17November2018/Mix_Fix_6h_CR-mtMask/", "cwFA_Mask", n = 500)
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



