library(SummarizedExperiment)
library(Matrix)
library(dplyr)
library(BuenColors)
library(ComplexHeatmap)
library(circlize)
library(data.table)
library(ggrepel)
source("../../global_functions/variant_calling.R")
load("../output/CD34_umap_embedding_granja.rda")

"%ni%" <- Negate("%in%")
perc.rank <- function(x) trunc(rank(x))/length(x)

# Import mgatk object
import_mgatk <- function(library){
  df <- fread(paste0("../output/barcode_qc/",library,".barcode_qc.tsv")) %>%
    filter(keep)
  se <- readRDS(paste0("../../../mtscATACpaper_large_data_files/source/mgatk_output/",library,"_v12-mtMask_mgatk.rds"))
  se_filt <- se[,as.character(df$sample)]
  colnames(se_filt) <- paste0(library, "_", colnames(se_filt))
  se_filt
}

#--------------------
# Call CD34 mutations
#--------------------
SE_CD34 <- cbind(import_mgatk("CD34_G10"), import_mgatk("CD34_H8"))
mut_se_CD34 <- call_mutations_mgatk(SE_CD34)
misc_df_CD34 <- data.frame(rowData(mut_se_CD34))

df <- misc_df_CD34 %>% filter(n_cells_conf_detected > 0 & n_cells_conf_detected <= 3 & n_cells_over_20 > 0) %>% filter(n_cells_over_5 <= 5) 
table(df$nucleotide)
dim(df)
rare_variant_cells <- names(cs)[cs > 0]

plot_df$rv_cell <- rownames(plot_df) %in% rare_variant_cells

plot_df %>% group_by(Clusters) %>% summarize(rv = sum(rv_cell), count = n()) %>%
  ungroup() %>% mutate(prop = rv/sum(rv), total = count/sum(count)) %>%
  mutate(FC = prop/total)


