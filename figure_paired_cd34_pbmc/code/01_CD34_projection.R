library(Matrix)
library(SummarizedExperiment)
library(tidyverse)
library(umap)
library(edgeR)
library(FNN)
library(matrixStats)
library(igraph)
library(BuenColors)
library(GenomicRanges)
set.seed(1)
"%ni%" <- Negate("%in%")

# Import granja data
if(FALSE){
  granja_mat <- fread("zcat < granja_cd34/GSE129785_scATAC-Hematopoiesis-CD34.mtx.gz", skip = 2)
  sm <- sparseMatrix(i = granja_mat[["V1"]], j = granja_mat[["V2"]], x = granja_mat[["V3"]])
  colData <- fread("granja_cd34/GSE129785_scATAC-Hematopoiesis-CD34.cell_barcodes.txt") %>% data.frame()
  peaks <- fread("granja_cd34/GSE129785_scATAC-Hematopoiesis-CD34.peaks.bed") %>%
    data.frame() %>% setnames(c('chr', 'start', 'end')) %>% makeGRangesFromDataFrame()
  SEall <- SummarizedExperiment(
    rowData = peaks,
    colData = colData, 
    assays = list(counts = sm)
  )
  cd34boo <- colData$Group %in% c("CD34_Progenitors_Rep1","CD34_Progenitors_Rep2")
  c1boo <-  colData$Group %in% c("BM_pDC", "CLP", "CMP", "GMP", "HSC", "LMPP", "MEP", "Monocytes", "MPP")
  SE_CD34 <- SEall[,cd34boo]
  SE_C1 <- SEall[,c1boo]
  saveRDS(SE_C1, file = "granja_cd34/granja_published_C1.rds")
  saveRDS(SE_CD34, file = "granja_cd34/granja_10X_CD34.rds")
  
}
SE <- readRDS("granja_cd34/granja_10X_CD34.rds")

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



#Set Clustering Parameters
nPCs1 <- 1:15
nPCs2 <- 1:15
resolution <- 0.8 #clustering resolution
nTop <- 10000 #number of variable peaks

#Run LSI 1st Iteration
lsi1 <- calcLSI(assay(SE), nComponents = 25, binarize = TRUE, nFeatures = NULL)
clust1 <- louvainIgraphClusters(lsi1[[1]], 10)

#Make Pseudo Bulk Library
message("Making PseudoBulk...")
clusterSums <- groupSums(mat = assay(SE), groups = clust1, sparse = TRUE) #Group Sums
logMat <- edgeR::cpm(clusterSums, log = TRUE, prior.count = 3) #log CPM matrix
varPeaks <- head(order(matrixStats::rowVars(logMat), decreasing = TRUE), nTop) #Top variable peaks

#Run LSI 2nd Iteration
lsi2 <- calcLSI(assay(SE)[varPeaks,,drop=FALSE], nComponents = 25, binarize = TRUE, nFeatures = NULL)
clust2 <- louvainIgraphClusters(lsi2[[1]], 10)
length(unique(clust2))

#UMAP
set.seed(1)
umap <- umap::umap(
  lsi2$matSVD[,nPCs1], 
  n_neighbors = 40, # original 55
  min_dist = 0.45, # original 0.45
  metric = "euclidean", 
  verbose = TRUE    )
set.seed(10)
plot_df <- data.frame(umap$layout, colData(SE), cl_clust = clust2)

p0 <- ggplot(plot_df, aes(x= X1, y = X2, color = Clusters)) +
   geom_point(size = 0.5) +
   labs(x = "UMAP1", y= "UMAP2", color = "") +
  pretty_plot() + L_border() + theme(legend.position = "bottom")

sel <- readRDS("granja_cd34/granja_published_C1.rds")
lsiProjection <- projectLSI((assay(sel)[varPeaks,]), lsi2)
umapProjection <- round(predict(umap, data.matrix(lsiProjection[,nPCs1])), 2)

projection_df <- data.frame(
  celltype = c(gsub("BM_", "", colData(sel)$Group), rep("none", dim(plot_df)[1])),
  umap1 = c(umapProjection[,1], plot_df$X1),
  umap2 = c(umapProjection[,2], plot_df$X2)
)

p1 <- ggplot(projection_df[dim(projection_df)[1]:1,], aes(x= umap1, y = umap2, color = celltype, label = celltype)) +
  geom_point(size = 0.5) +
  labs(x = "UMAP1", y= "UMAP2", color = "C1 FACS ") +
  pretty_plot() + L_border() + theme(legend.position = "bottom") +
  scale_color_manual(values = c(ejc_color_maps, "none" = "lightgrey", "Monocytes" = "orange2"))
p1

# Now project the Pearson cells
pearson_cd34 <- readRDS("rds_SE/Pearson_CD34_QC.rds")
lsiProjectionPearson <- projectLSI((assay(pearson_cd34)[varPeaks,]), lsi2)
umapProjectionPearson <- round(predict(umap, data.matrix(lsiProjectionPearson[,1:25])), 2)

projection_df_p <- data.frame(
  celltype = c(rep("Pearson", length(umapProjectionPearson[,1])), rep("none", dim(plot_df)[1])),
  umap1 = c(umapProjectionPearson[,1], plot_df$X1),
  umap2 = c(umapProjectionPearson[,2], plot_df$X2),
  bc= c(rownames(umapProjectionPearson), plot_df$Group_Barcode)
)

het_df <- fread("data/PT3_all_libraries.tsv") %>% filter(cell_clust == "BM_CD34")
mds_del <- fread("../del7_CD34.tsv")

mdf <- left_join(projection_df_p, het_df, by = c("bc" = "cell_id"))
mdf2 <- left_join(mdf, mds_del, by = c("bc" = "barcodes"))
mdf2$deleted <- ifelse(mdf2$percent_del < 25, "Monosomy7", "WT7")

p2 <- ggplot(mdf[dim(mdf)[1]:1,], aes(x= umap1, y = umap2, color = celltype, label = celltype)) +
  geom_point(size = 0.5) +
  labs(x = "UMAP1", y= "UMAP2", color = "Sample") +
  pretty_plot() + L_border() + theme(legend.position = "bottom") +
  scale_color_manual(values = c("lightgrey", "firebrick"))

p3 <- ggplot(mdf[dim(mdf)[1]:1,], aes(x= umap1, y = umap2, color = corr_heteroplasmy, label = corr_heteroplasmy)) +
  geom_point(size = 0.5) +
  scale_color_gradientn(colors = jdb_palette("solar_extra"), na.value = "lightgrey") +
  labs(x = "UMAP1", y= "UMAP2", color = "Deletion %") +
  pretty_plot() + L_border() + theme(legend.position = "bottom")

p4 <- ggplot(mdf2[dim(mdf2)[1]:1,], aes(x= umap1, y = umap2, color = deleted)) +
  geom_point(size = 0.5) +
  labs(x = "UMAP1", y= "UMAP2", color = "") +
  pretty_plot() + L_border() + theme(legend.position = "bottom") +
  scale_color_manual(values = c("firebrick", "dodgerblue3"), na.value = "lightgrey")

ggplot(mdf2 %>% filter(celltype == "Pearson"), aes(x = deleted, y = corr_heteroplasmy)) +
  geom_beeswarm(size = 0.2) + geom_boxplot(color = "firebrick", fill = NA) +
  pretty_plot() + L_border() + theme(legend.position = "bottom") +
  labs(x = "", y = "Heteroplasmy")

cowplot::ggsave(cowplot::plot_grid(p0, p1, p2, p3, p4, p5, nrow = 3), 
                file = "Pearson_CD34_embedding.pdf", width = 8, height = 13)

