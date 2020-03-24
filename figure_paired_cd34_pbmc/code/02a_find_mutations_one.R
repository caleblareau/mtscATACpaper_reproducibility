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

df <- misc_df_CD34 %>% filter(n_cells_conf_detected == 1 & n_cells_over_20 > 0) %>% filter(n_cells_over_5 <= 3) 
table(df$nucleotide)
