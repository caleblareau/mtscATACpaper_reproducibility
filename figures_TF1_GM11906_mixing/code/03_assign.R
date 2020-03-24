library(SummarizedExperiment)
library(Matrix)
library(dplyr)
library(BuenColors)
library(ComplexHeatmap)
library(circlize)
"%ni%" <- Negate("%in%")
source("02a_assign_pull_function.R")

df1 <- assign_pull(readRDS("../../../mtscATACpaper_large_data_files/source/mgatk_output/Mix_Fix_1h_CR-mtMask_mgatk.rds"))
df6 <- assign_pull(readRDS("../../../mtscATACpaper_large_data_files/source/mgatk_output/Mix_Fix_6h_CR-mtMask_mgatk.rds"))
dfMix1 <- assign_pull(readRDS("../../../mtscATACpaper_large_data_files/source/mgatk_output/Mix_1_CR-mtMask_mgatk.rds"))

write.table(df1, file = "../output/data1_meta.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(df6, file = "../output/data6_meta.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(dfMix1, file = "../output/data_mix1_meta.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
