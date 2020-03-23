#Clustering and scATAC-seq UMAP for Hematopoiesis data
#06/02/19
#Cite Granja*, Klemm*, Mcginnis* et al. 
#A single cell framework for multi-omic analysis of disease identifies 
#malignant regulatory signatures in mixed phenotype acute leukemia (2019)
#Created by Jeffrey Granja
library(Matrix)
library(SummarizedExperiment)
library(tidyverse)
library(umap)
library(edgeR)
library(FNN)
library(matrixStats)
library(igraph)
library(BuenColors)
set.seed(1)

####################################################
#Functions
####################################################

#Binarize Sparse Matrix
binarizeMat <- function(mat){
  mat@x[mat@x > 0] <- 1
  mat
}

#LSI Adapted from fly-atac with information for re-projection analyses
calcLSI <- function(mat, nComponents = 50, binarize = TRUE, nFeatures = NULL, log = TRUE){
  
  set.seed(1)
  
  #TF IDF LSI adapted from flyATAC
  if(binarize){
    message(paste0("Binarizing matrix..."))
    mat@x[mat@x > 0] <- 1 
  }
  
  if(!is.null(nFeatures)){
    message(paste0("Getting top ", nFeatures, " features..."))
    idx <- head(order(Matrix::rowSums(mat), decreasing = TRUE), nFeatures)
    mat <- mat[idx,] 
  }else{
    idx <- which(Matrix::rowSums(mat) > 0)
    mat <- mat[idx,]
  }
  
  #Calc RowSums and ColSums
  colSm <- Matrix::colSums(mat)
  rowSm <- Matrix::rowSums(mat)
  
  #Calc TF IDF
  message("Computing Term Frequency IDF...")
  freqs <- t(t(mat)/colSm)
  idf   <- as(log(1 + ncol(mat) / rowSm), "sparseVector")
  if(log){
    freqs@x = log1p(freqs@x * 100000) # equivalent to adding a small pseudocount, but without making matrix dense
    tfidf <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% freqs
    
  } else{
    tfidf <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% freqs
  }
  
  
  #Calc SVD then LSI
  message("Computing SVD using irlba...")
  svd <- irlba::irlba(tfidf, nComponents, nComponents)
  svdDiag <- matrix(0, nrow=nComponents, ncol=nComponents)
  diag(svdDiag) <- svd$d
  matSVD <- t(svdDiag %*% t(svd$v))
  rownames(matSVD) <- colnames(mat)
  colnames(matSVD) <- paste0("PC",seq_len(ncol(matSVD)))
  
  #Return Object
  out <- list(
    matSVD = matSVD, 
    rowSm = rowSm, 
    colSm = colSm, 
    idx = idx, 
    svd = svd, 
    binarize = binarize, 
    nComponents = nComponents,
    date = Sys.Date(),
    seed = 1)
  
  out
  
}

projectLSI <- function(mat, lsi, log = TRUE){   
  
  #Get Same Features
  mat <- mat[lsi$idx,]
  if(lsi$binarize){
    message(paste0("Binarizing matrix..."))
    mat@x[mat@x > 0] <- 1       
  }
  
  #Calc TF IDF
  rowsToZero <- which(lsi$rowSm == 0)
  setToZero <- which((mat@i + 1) %in% rowsToZero)
  if(length(setToZero) > 0){
    mat@x[setToZero] <- 0
  }
  
  message("Computing Term Frequency IDF...")
  freqs <- t(t(mat)/Matrix::colSums(mat))
  idf   <- as(log(1 + length(lsi$colSm) / lsi$rowSm), "sparseVector")
  
  if(log){
    freqs@x = log1p(freqs@x * 100000) # equivalent to adding a small pseudocount, but without making matrix dense
    tfidf <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% freqs
    
  } else{
    tfidf <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% freqs
  }
  
  if(length(Matrix::which(is.na(tfidf),arr.ind=TRUE)) > 0){
    tfidf[Matrix::which(is.na(tfidf),arr.ind=TRUE)] <- 0 #weird Inf * 0
  }
  
  #Calc V
  V <- t(tfidf) %*% lsi$svd$u %*% diag(1/lsi$svd$d)
  
  #Calc SVD then LSI
  message("Computing SVD using irlba...")
  svdDiag <- matrix(0, nrow=lsi$nComponents, ncol=lsi$nComponents)
  diag(svdDiag) <- lsi$svd$d
  matSVD <- t(svdDiag %*% t(V))
  rownames(matSVD) <- colnames(mat)
  colnames(matSVD) <- paste0("PC",seq_len(ncol(matSVD)))
  
  return(matSVD)
  
}

#Helper function for summing sparse matrix groups
groupSums <- function (mat, groups = NULL, na.rm = TRUE, sparse = FALSE){
  stopifnot(!is.null(groups))
  stopifnot(length(groups) == ncol(mat))
  gm <- lapply(unique(groups), function(x) {
    if (sparse) {
      Matrix::rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
    }
    else {
      rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
    }
  }) %>% Reduce("cbind", .)
  colnames(gm) <- unique(groups)
  return(gm)
}

#Seurat SNN
seuratSNN <- function(matSVD, dims.use = 1:50, print.output = TRUE, ...){
  set.seed(1)
  message("Making Seurat Object...")
  mat <- matrix(rnorm(nrow(matSVD) * 3, 1000), ncol = nrow(matSVD), nrow = 3)
  colnames(mat) <- rownames(matSVD)
  obj <- Seurat::CreateSeuratObject(mat, project='scATAC', min.cells=0, min.genes=0)
  obj <- Seurat::SetDimReduction(object = obj, reduction.type = "pca", slot = "cell.embeddings", new.data = matSVD)
  obj <- Seurat::SetDimReduction(object = obj, reduction.type = "pca", slot = "key", new.data = "PC")
  obj <- Seurat::FindClusters(object = obj, reduction.type = "pca", dims.use = dims.use, print.output = print.output, ...)
  clust <- obj@meta.data[,ncol(obj@meta.data)]
  paste0("mc",match(clust, unique(clust)))
}

louvainIgraphClusters <- function(matSVD, k){
  set.seed(1)
  knn <- FNN::get.knn(matSVD, algo="kd_tree", k = k)[["nn.index"]]
  igraphObj <- igraph::graph_from_adjacency_matrix(igraph::get.adjacency(igraph::graph.edgelist(data.matrix(reshape2::melt(knn)[,c("Var1", "value")]), directed=FALSE)), mode = "undirected")
  clusters <- membership(cluster_louvain(igraphObj))
  return(paste0("mc", as.character(clusters)))
}

umapIgraphClusters <- function(umap, k){
  set.seed(1)
  knn <- FNN::get.knn(umap, algo="kd_tree", k = k)[["nn.index"]]
  igraphObj <- igraph::graph_from_adjacency_matrix(igraph::get.adjacency(igraph::graph.edgelist(data.matrix(reshape2::melt(knn)[,c("Var1", "value")]), directed=FALSE)), mode = "undirected")
  clusters <- membership(cluster_louvain(igraphObj))
  return(paste0("mc", as.character(clusters)))
}


build_means_for_clusters <- function(LSI_mat, clusters){
  sapply(sort(unique(clusters)), function(cluster_i){
    colMeans(LSI_mat[clusters == cluster_i,])
  })
}

# A: A matrix of clusters x PCs
# B: A matrix of cells x PCs
classify_from_reference <- function(A, B){
  
  # Get the pair-wise euclidean distance msot-alike
  euklDist <- t(sqrt(apply(array(apply(B,1,function(x){(x-t(A))^2}),c(ncol(A),nrow(A),nrow(B))),2:3,sum)))
  colnames(euklDist) <- rownames(A)
  colnames(euklDist)[max.col(-1*euklDist, 'first')]  -> vec
  return(vec)
}

####################################################
#Input Data
####################################################
#Read in Summarized Experiment
#Please Note Code here has been modified to work with finalized summarized experiment

#Reference Summarized Experiment
#Contains Peaks for Reference Hematopoiesis only
seReference <-readRDS("rds_SE/PBMC_5k_nextgem_QC.rds")

#Set Clustering Parameters
nPCs1 <- 1:25
nPCs2 <- 1:25
resolution <- 0.8 #clustering resolution
nTop <- 25000 #number of variable peaks

#Create Matrix
mat <- assay(seReference)

#Run LSI 1st Iteration
lsi1 <- calcLSI(mat, nComponents = 50, binarize = TRUE, nFeatures = NULL)
clust1 <- louvainIgraphClusters(lsi1[[1]], 10)

#Make Pseudo Bulk Library
message("Making PseudoBulk...")
mat <- mat[,rownames(lsi1[[1]]), drop = FALSE] #sometimes cells are filtered
mat@x[mat@x > 0] <- 1 #binarize
clusterSums <- groupSums(mat = mat, groups = clust1, sparse = TRUE) #Group Sums
logMat <- edgeR::cpm(clusterSums, log = TRUE, prior.count = 3) #log CPM matrix
varPeaks <- head(order(matrixStats::rowVars(logMat), decreasing = TRUE), nTop) #Top variable peaks

#Run LSI 2nd Iteration
lsi2 <- calcLSI(mat[varPeaks,,drop=FALSE], nComponents = 50, binarize = TRUE, nFeatures = NULL)
clust2 <- louvainIgraphClusters(lsi2[[1]], 20)
length(unique(clust2))
means_lsi <- build_means_for_clusters(lsi2$matSVD, clust2)
lsiPeaks <- varPeaks
lsiReference <- lsi2

#UMAP
set.seed(1)
umap <- umap::umap(
  lsi2$matSVD[,1:25], 
  n_neighbors = 30, # original 55
  min_dist = 0.4, # original 0.45
  metric = "cosine", 
  verbose = TRUE    )

clust2 <- umapIgraphClusters(umap$layout, k = 100)
length(unique(clust2))
means_lsi <- build_means_for_clusters(lsi2$matSVD, clust2)

ggplot(data.frame(umap$layout, clusters = clust2), aes(x = X1, y = X2, color = clust2)) + geom_point() +
  scale_color_manual(values = (jdb_palette("lawhoops")))


# Make classification data.frame
important_df <- data.frame(
  barcode = c(colnames(seReference)),
  batch = c(rep("base", dim(seReference)[2])),
  umap1 = c(umap$layout[,1]),
  umap2 = c(umap$layout[,2]),
  cluster = c(clust2)
)

annovec <- c("myeloid", "myeloid", "CD4", "CD4", "Bcell", "CD8", "myeloid", "myeloid", "myeloid", "NKcell", "CD4", "CD8", "Bcell")
names(annovec) <- paste0("mc", as.character(1:13))
important_df$cell_cluster <- annovec[as.character(important_df$cluster)]
important_df %>% group_by(cell_cluster) %>%
  summarize(umap1= mean(umap1), umap2 = mean(umap2)) %>% ggplot(aes(x = umap1, y = umap2, label = cell_cluster)) +
  geom_text()

save(important_df, umap, lsiPeaks, means_lsi, lsiReference, file = "15November2019_projection_5k_CLL.rda")

p1 <- ggplot(important_df %>% filter(batch == "base"), aes(x = umap1, y = umap2, color = cell_cluster)) +
  geom_point_rast(size = 0.05, raster.dpi = 1000) +
  pretty_plot(fontsize = 8) + L_border() + labs(x = "UMAP1", y = "UMAP2") +
  theme(legend.position = "none") + scale_color_manual(values = c("Bcell" = "#BA7FD0", "CD4" ="#0081C9", "CD8"="#001588", "NKcell"="#490C65", "myeloid"="#FF5A00"))


cowplot::ggsave(p1, file = "plots/base_umap_embedding.pdf", width = 1.7, height = 1.7)

if(FALSE){
  frags_base <- fread("public_data/atac_pbmc_5k_nextgem_fragments.tsv.gz")

  clusterss <- sort(unique(as.character(important_df$cluster)))
  lapply(clusterss, function(cc){
    
    bc_base <- important_df %>% filter(cluster == cc & batch == "base") %>% pull(barcode)

    # Setup up low heteroplasmy
    gr <- frags_base[V4 %in% bc_base]  %>% data.frame() %>%
      setnames(c("chr", "start", "end", "bc", "pcrn")) %>% 
      makeGRangesFromDataFrame()
    reads_coverage <- coverage(gr)/length(gr)*1000000
    
    rtracklayer::export.bw(reads_coverage, con = paste0("bw/", "5k_",cc,".bw"))
    cc
  })
}

