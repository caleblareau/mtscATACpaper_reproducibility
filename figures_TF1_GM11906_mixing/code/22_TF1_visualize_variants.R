library(SummarizedExperiment)
library(Matrix)
library(dplyr)
library(BuenColors)
library(ComplexHeatmap)
library(circlize)
library(irlba)
library(Seurat)

"%ni%" <- Negate("%in%")
perc.rank <- function(x) trunc(rank(x))/length(x)
source("../../global_functions/get_allele_freq_mat.R")

#Import data
# Import mgatk object
SE <- cbind(readRDS("../../../mtscATACpaper_large_data_files/source/mgatk_output/Mix_Fix_1h_CR-mtMask_mgatk.rds"), 
            readRDS("../../../mtscATACpaper_large_data_files/source/mgatk_output/Mix_Fix_6h_CR-mtMask_mgatk.rds"))

# Find TF1 cells
rbind(read.table("../output/data1_meta.tsv", header = TRUE), 
      read.table("../output/data6_meta.tsv", header = TRUE)) %>% filter(assign == "TF1") %>%
  filter(mean_cov > 50) -> TF_cells_df

SE <- SE[,as.character(TF_cells_df$cell_id)]
mmat <- computeAFMutMatrix(SE)
  
# Get mutation matrix and threshold
vars <- read.table("../output/TF1_VMR_processed_variants_ndetect.tsv", header = FALSE, stringsAsFactors = FALSE)[,1]
afp <- mmat[vars,]
#afp[afp>0.2] <- 0.2

# Get clusters
seuratSNN <- function(matSVD, resolution = 0.8, k.param = 20){ # consistent with seurat defaults
  set.seed(1)
  rownames(matSVD) <- make.unique(rownames(matSVD))
  obj <- FindNeighbors(matSVD, k.param = k.param, annoy.metric = "cosine")
  clusters <- FindClusters(object = obj$snn, resolution = resolution)
  return(as.character(clusters[,1]))
}

cl <- seuratSNN(sqrt(t(afp)))
table(cl)

afp_svd <- prcomp_irlba(sqrt(afp), n = 10)
mat <- afp_svd$rotation
rownames(mat) <- colnames(afp)

clusters <- seuratSNN(mat)
table(clusters)

set.seed(2)
names_clusters <- unique((cl))
vec_go <- sample(as.character(jdb_palettes[["lawhoops"]]))[1:length(names_clusters)]
names(vec_go) <- names_clusters

# Make data.frame for stuff
df <- data.frame(
  cell_id = colnames(afp), 
  cluster_id = as.character(cl)
) %>% arrange(desc(cl))

ha_col <- HeatmapAnnotation(cell = as.character(df$cluster_id),
                            col = list(cell = vec_go))

# Manually configure variant order
var_order <- c("8002C>T", "7789G>C","4037G>A","14569G>A","2954C>A", "5970G>A",'7521G>A',
               "11196G>A", "309C>T",  "7847G>A", 
               "7256C>T",  "7598G>A",
               "13708G>A", "627G>A", 
               "14905G>A", "5560G>A", "1793G>A", "9493G>A","12433C>T","12316G>A",
               "3901G>A", "10552C>T",
               "7316G>A", "1888G>A", "6021G>A", "6713C>T",  "2643G>A",
               "709G>A", "7276C>T",
               "5790C>T", "2002G>A", "5545C>A", "12012C>T",
               "8684C>T", "10458C>T", "15765G>A", "2573G>A","6027G>A","2040G>A","1045G>A",
               "5324C>T","2789C>T",
               "7457G>A", "10068G>A", "6762G>A", "12762C>A", "14559C>T","15646C>T"
)
var_order[duplicated(var_order)]
setdiff(vars, var_order)
setdiff(var_order, vars)

aftree <- afp
afp[afp < 0.01] <- 0
afp[afp > 0.1] <- 0.1

pdf(paste0("../plots/final_TF1.pdf"), width=4, height=3)
hm <- Heatmap((data.matrix(afp)[var_order,as.character(df$cell_id)]),  # var_order
              col=as.character(jdb_palette("solar_rojos",type="continuous")),
              show_row_names = TRUE, 
              top_annotation=ha_col,
              cluster_columns = FALSE,
              name = "AF",
              row_names_gp = gpar(fontsize = 4),
              cluster_rows = FALSE, 
              show_column_names = FALSE)
hm 
dev.off()


png(paste0("../plots/final_TF1_nosqrt_10.png"), width=4, height=3, unit = "in", res = 500)
hm <- Heatmap((data.matrix(afp)[var_order,as.character(df$cell_id)]),
              col=as.character(jdb_palette("solar_rojos",type="continuous")),
              show_row_names = TRUE, 
              top_annotation=ha_col,
              cluster_columns = FALSE,
              name = "AF",
              row_names_gp = gpar(fontsize = 4),
              cluster_rows = FALSE, 
              show_column_names = FALSE)
hm 
dev.off()


library(ggplot2)
library(phangorn)
library(ape)
library(ggtree)

# Get group means
matty <- sapply(names(vec_go), function(cluster){
  cells <- df %>% filter(cluster_id == cluster) %>% pull(cell_id) %>% as.character()
  Matrix::rowMeans(sqrt(afp[,cells]))
})

mito.nj <- phangorn::NJ(dist(t(sqrt(matty))))
pdf("../plots/tree_me.pdf", width = 4, height = 4)
plot(mito.nj)
dev.off()
