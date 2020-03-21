library(SummarizedExperiment)
library(Matrix)
library(dplyr)
library(BuenColors)
library(ComplexHeatmap)
library(circlize)
"%ni%" <- Negate("%in%")
source("../../global_functions/get_allele_freq_mat.R")

# This is a first pass QC script to find homoplasmic variants present in either cell line
# Not used in the final figures for the manuscript

# Import AFs
SE <- readRDS("../../../mtscATACpaper_large_data_files/source/mgatk_output/Mix_Fix_6h_CR-mtMask_mgatk.rds")
filt <- colData(SE)$depth >= 50
SE2 <- SE[,filt]
af <- computeAFMutMatrix( SE2 )

x <- (rowSums(af > 0.99) > 5) & (rowSums(af < 0.005) > 5)
sum(x)
xx <- (x  & (rownames(af) %ni% c("3107N>C", "3107N>A", "185G>A", "3107N>T", "13708G>A", "8860A>G") )) | rownames(af) %in% "8344A>G"
afp <- af[xx,]

perc.rank <- function(x) trunc(rank(x))/length(x)
coverage_quantile <- perc.rank(colData(SE2)$depth)

ha_col <- HeatmapAnnotation(df = data.frame(CovRank = coverage_quantile),
                            col = list(CovRank = colorRamp2(c(0,1), c("white", "dodgerblue3"))))

pdf(paste0("../plots/quick_mix_fix_6_find_homoplasmic.pdf"), width=10, height=10)
hm <- Heatmap(data.matrix(afp), 
              col=as.character(jdb_palette("solar_rojos",type="continuous")),
              show_row_names = TRUE, 
              cluster_columns = TRUE,
              top_annotation=ha_col,
              name = "AF",
              cluster_rows = TRUE, 
              show_column_names = FALSE)
hm
dev.off()
dim(af)

