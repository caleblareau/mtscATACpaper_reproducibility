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

# Import signac procssed stuff
d <- readRDS("../output/21March2020_signac_process.rds")
d$barcode <- rownames(d)

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

# Cap for visualiation
zscore[zscore > 3] <- 3
zscore[zscore < -3] <- -3

keep <- TRUE
ordering <- factor(region_meta[keep,1], levels = unique(region_meta[keep,1]))

colVec <- c("dodgerblue3", "purple4")
names(colVec) <- as.character(c(0,1))

d_ss <- d %>% filter(seurat_clusters %in% c(0,1)) %>% arrange(seurat_clusters)
ha1 <- HeatmapAnnotation(df = data.frame(clusters = as.character(d_ss$seurat_clusters)), 
                         col = list(clusters = colVec)
)

pdf(paste0("../plots/CRC_tumor_CNV.pdf"), width=3, height=5)
Heatmap(zscore[keep,as.character(d_ss$barcode)], 
        col=as.character(jdb_palette("solar_basic",type="continuous")),
        show_row_names = FALSE, 
        cluster_columns = FALSE,
        name = "CNV zscore",
        row_names_gp = gpar(fontsize = 3),
        cluster_rows = FALSE, 
        split = ordering,
        top_annotation = ha1,
        show_column_names = FALSE)
dev.off()



ggplot(d, aes(y = pct_reads_in_DNase, x =log10_mtDNA_depth, color = seurat_clusters)) +
  geom_point() + scale_color_manual(values = c("dodgerblue3", "purple4", "forestgreen", "orange3", "firebrick", "orange", "green2")) 

ggplot(d, aes(x = pct_reads_in_DNase, y = chr12, color  = seurat_clusters)) + geom_point() +
  scale_color_manual(values = c("dodgerblue3", "purple4", "forestgreen", "orange3", "firebrick", "orange", "green2")) 

