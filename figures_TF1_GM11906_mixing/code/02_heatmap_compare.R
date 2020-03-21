library(SummarizedExperiment)
library(Matrix)
library(dplyr)
library(BuenColors)
library(ComplexHeatmap)
library(circlize)
"%ni%" <- Negate("%in%")

source("../../global_functions/get_allele_freq_mat.R")
source("02a_assign_pull_function.R")

# Import AFs
SE <- readRDS("../data/big_gi/Mix_Fix_1h_CR-mtMask_mgatk.rds")
filt <- colData(SE)$depth > 20
SE2 <- SE[,filt]
af <- computeAFMutMatrix( SE2 )
assign_df <- assign_pull(SE2) %>% arrange(tf1_mean)
muts <- rev(read.table("../data/variant_order.txt", header = FALSE)[,1] %>% as.character())
muts <- muts[muts %ni% c("12705C>T")]
afp <- af[muts,as.character(assign_df$cell_id)]
png(paste0("../output/mix_good.png"), width=10, height=10, units = "in", res = 300)
hm <- Heatmap(data.matrix(afp), 
              col=as.character(jdb_palette("solar_rojos",type="continuous")),
              show_row_names = FALSE, 
              cluster_columns = FALSE,
              name = "AF",
              cluster_rows = FALSE, 
              show_column_names = FALSE)
hm
dev.off()

