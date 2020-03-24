library(SummarizedExperiment)
library(Matrix)
library(dplyr)
library(BuenColors)
library(ComplexHeatmap)
library(circlize)
"%ni%" <- Negate("%in%")

source("../../global_functions/get_allele_freq_mat.R")
source("02a_assign_pull_function.R")

# Import data and compute allele frequencies
SE <- readRDS("../../../mtscATACpaper_large_data_files/source/mgatk_output/Mix_Fix_1h_CR-mtMask_mgatk.rds")
filt <- colData(SE)$depth > 20
SE2 <- SE[,filt]
af <- computeAFMutMatrix( SE2 )
assign_df <- assign_pull(SE2) %>% arrange(tf1_mean)
muts <- rev(read.table("../data/variant_order.txt", header = FALSE)[,1] %>% as.character())
muts <- muts[muts %ni% c("12705C>T")]
afp <- af[muts,as.character(assign_df$cell_id)]

# Export heatmap
png(paste0("../plots/mix_good_heatmap.png"), width=10, height=10, units = "in", res = 300)
hm <- Heatmap(data.matrix(afp), 
              col=as.character(jdb_palette("solar_rojos",type="continuous")),
              show_row_names = FALSE, 
              cluster_columns = FALSE,
              name = "AF",
              cluster_rows = FALSE, 
              show_column_names = FALSE)
hm
dev.off()

# Do the same for the no fixation sample
SE <- readRDS("../../../mtscATACpaper_large_data_files/source/mgatk_output/Mix_2_NP40_CR-mtMask_mgatk.rds")
filt <- colData(SE)$depth > 20
SE2 <- SE[,filt]
af <- computeAFMutMatrix( SE2 )
assign_df <- assign_pull(SE2) %>% arrange(tf1_mean)
muts <- rev(read.table("../data/variant_order.txt", header = FALSE)[,1] %>% as.character())
muts <- muts[muts %ni% c("12705C>T")]
afp <- af[muts,as.character(assign_df$cell_id)]

# Export heatmap
png(paste0("../plots/mix_bad_heatmap.png"), width=10, height=10, units = "in", res = 300)
hm <- Heatmap(data.matrix(afp), 
              col=as.character(jdb_palette("solar_rojos",type="continuous")),
              show_row_names = FALSE, 
              cluster_columns = FALSE,
              name = "AF",
              cluster_rows = FALSE, 
              show_column_names = FALSE)
hm
dev.off()


