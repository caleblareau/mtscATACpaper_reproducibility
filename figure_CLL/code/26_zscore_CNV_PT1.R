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


# Get the typical profile of a fragment distribution
base_pbmcs <- readRDS("../../cnv_compute/output/atac_public_pbmcs_cnv.rds"); base_pbmcs[is.na(base_pbmcs)] <- 0

cpm_norm <- (t(t(base_pbmcs)/colSums(base_pbmcs)) * 100000)
row_means <- rowMeans(cpm_norm)
row_std <- sqrt(rowVars(cpm_norm))

mat <- readRDS("../../cnv_compute/output/CLL_PT1_CD19pos_cnv.rds"); mat[is.na(mat)] <- 0
clone_definition <- readRDS("../output/PT1_clone_definition.rds")
bcs <- as.character(clone_definition$cell_id)
mat <- mat[,bcs]
makeZscoreMat <- function(mat){
  mat_cpm_norm <- (t(t(mat)/colSums(mat)) * 100000)
  zscore_mat <- (mat_cpm_norm - row_means)/row_std
}

makeLog2Mat <- function(mat){
  mat_cpm_norm <- (t(t(mat)/colSums(mat)) * 100000)
  zscore_mat <- log2(mat_cpm_norm/(row_medians + 1))
  zscore_mat
}

zscore <- makeZscoreMat(mat)
score_bs <- makeZscoreMat(base_pbmcs)
region_meta <- str_split_fixed(rownames(zscore), "_", 4)

clone_definition$zscore <- colMeans(zscore[region_meta[,1] == "chr12",])
baseline_df <- data.frame(
  cell_id = "whatever",
  cluster_id = c("b"), 
  zscore = colMeans(score_bs[region_meta[,1] == "chr12",])
)

vec_go <- c("firebrick", jdb_palettes[["flame_light"]],  "grey")

p1 <- ggplot(rbind(clone_definition, baseline_df), aes(x = 1, y = zscore, color = cluster_id)) +
  geom_boxplot(outlier.shape = NA, width = 0.1, fill = NA) +
  scale_y_continuous(limits = c(-1, 3)) + labs(x = "Cell clone", y = "chr12 CNV z-score") +
  pretty_plot() + L_border() + scale_color_manual(values = vec_go) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(legend.position = "none")
  
cowplot::ggsave2(p1, file = "../plots/PT1_CNVplot.pdf", width = 2, height = 1.7)

# Compute fraction of chromsome 12 reads

# normal
f12n <- colSums(mat[region_meta[,1] == "chr12",])/colSums(mat)
f12c <- colSums(base_pbmcs[region_meta[,1] == "chr12",])/colSums(base_pbmcs)

mean(f12n)/mean(f12c)
