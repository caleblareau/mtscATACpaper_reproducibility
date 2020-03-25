library(SummarizedExperiment)
library(Matrix)
library(dplyr)
library(BuenColors)
library(ComplexHeatmap)
library(circlize)
library(irlba)
library(Seurat)
library(stringr)

"%ni%" <- Negate("%in%")
perc.rank <- function(x) trunc(rank(x))/length(x)
source("../../global_functions/get_allele_freq_mat.R")

#Import data
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

# Get clusters
seuratSNN <- function(matSVD, resolution = 1, k.param = 10){ 
  set.seed(1)
  rownames(matSVD) <- make.unique(rownames(matSVD))
  obj <- FindNeighbors(matSVD, k.param = k.param, annoy.metric = "cosine")
  clusters <- FindClusters(object = obj$snn, resolution = resolution)
  return(as.character(clusters[,1]))
}

cl <- seuratSNN(sqrt(t(afp)))
cl <- str_pad(cl, 2, pad = "0")

set.seed(2)
names_clusters <- unique(cl)
colors <- c("#1f77b4","#d62728","#2ca02c","#ff7f0e","#9467bd","#8c564b","#e377c2","#7f7f7f",
            "#bcbd22","#17becf","#ad494a","#e7ba52","#8ca252","#756bb1","#636363","#aec7e8")
vec_go <- colors[1:length(names_clusters)]
names(vec_go) <- sort(names_clusters)

# Make data.frame for stuff
df <- data.frame(
  cell_id = colnames(afp), 
  cluster_id = as.character(cl)
) %>% arrange(cl)

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


library(phangorn)
library(ape)

# Get group means 
matty <- sapply(names(vec_go), function(cluster){
  cells <- df %>% filter(cluster_id == cluster) %>% pull(cell_id) %>% as.character()
  Matrix::rowMeans(sqrt(afp[,cells]))
})

# Do cosine distance; note that we used sqrt transformation already 
# when creating the pseudo bulk-cell matrix
mito.hc <- hclust(dist(lsa::cosine((matty))))
plot(mito.hc)

pdf("../plots/hierarchical_tree_TF1.pdf", width = 5, height = 5)
plot(mito.hc)
dev.off()
