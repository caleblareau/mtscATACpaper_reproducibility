library(SummarizedExperiment)
library(BuenColors)
library(dplyr)
library(Matrix)
library(BuenColors)
library(data.table)
library(ggbeeswarm)
library(dplyr)
library(matrixStats)
library(ComplexHeatmap)
library(circlize)
library(stringr)

# Compute baseline features per bin
base_pbmcs <- readRDS("../../cnv_compute/output/atac_public_pbmcs_cnv.rds"); base_pbmcs[is.na(base_pbmcs)] <- 0
cpm_norm <- (t(t(base_pbmcs)/colSums(base_pbmcs)) * 100000)
row_means <- rowMeans(cpm_norm)
row_std <- sqrt(rowVars(cpm_norm))

# Subset observed new matrix to valid barcodes 
mat <- readRDS("../../cnv_compute/output/CRC_tumor_cnv.rds"); mat[is.na(mat)] <- 0
bcs <- fread("../output/CRC_filtered_barcodes.tsv", header = FALSE)[[1]]
mat <- mat[,bcs]

# Functions to make Z score and log2 change w.r.t. baseline from 10x
makeZscoreMat <- function(mat){
  mat_cpm_norm <- (t(t(mat)/colSums(mat)) * 100000)
  zscore_mat <- (mat_cpm_norm - row_means)/row_std
}

makeLog2Mat <- function(mat){
  mat_cpm_norm <- (t(t(mat)/colSums(mat)) * 100000)
  zscore_mat <- log2(mat_cpm_norm/(row_medians + 1))
  zscore_mat
}

# Compute Zscores
zscore <- makeZscoreMat(mat)
score_base <- makeZscoreMat(base_pbmcs)
region_meta <- str_split_fixed(rownames(zscore), "_", 4)


ggplot(d, aes(x = pct_reads_in_DNase, y = chr12, color  = seurat_clusters)) + geom_point() +
  scale_color_manual(values = c("dodgerblue3", "purple4", "forestgreen", "orange3", "firebrick", "orange", "green2")) 


zscore[zscore > 3] <- 3
zscore[zscore < -3] <- -3

keep <- TRUE
ordering <- factor(region_meta[keep,1], levels = unique(region_meta[keep,1]))

colVec <- c("dodgerblue3", "purple4", "forestgreen", "orange3", "firebrick", "orange", "green2")
names(colVec) <- as.character(0:6)

ha1 <- HeatmapAnnotation(df = data.frame(clusters = as.character(mito_df$seurat_clusters)), 
                         col = list(clusters = colVec)
)

png(paste0("../plots/CRC_tumor_CNV.png"), width=5, height=10, units = "in", res = 500)
Heatmap(zscore[keep,], 
        col=as.character(jdb_palette("solar_basic",type="continuous")),
        show_row_names = FALSE, 
        cluster_columns = TRUE,
        name = "CNV zscore",
        row_names_gp = gpar(fontsize = 4),
        cluster_rows = FALSE, 
        split = ordering,
        top_annotation = ha1,
        show_column_names = FALSE)
dev.off()

ggplot(d, aes(y = pct_reads_in_DNase, x =log10_mtDNA_depth, color = seurat_clusters)) +
  geom_point() + scale_color_manual(values = c("dodgerblue3", "purple4", "forestgreen", "orange3", "firebrick", "orange", "green2")) 
