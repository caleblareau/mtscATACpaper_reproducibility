library(Matrix)
library(SummarizedExperiment)
library(tidyverse)
library(umap)
library(edgeR)
library(FNN)
library(matrixStats)
library(igraph)
library(BuenColors)
library(ggrastr)
set.seed(1)
"%ni%" <- Negate("%in%")

# Import the peaks x cells summarized experiment object
import_SE <- function(library){
  x <- readRDS(paste0("../../../mtscATACpaper_large_data_files/intermediate/",library,".rds"))
  colData(x)$library <- library
  colnames(x) <- paste0(library,"-", colnames(x))
  return(x)
}

# Import deviations from chromVAR for annotation
import_CV <- function(library){
  x <- readRDS(paste0("../output/chromVAR/",library,"_tf_deviations.rds"))
  colData(x)$library <- library
  colnames(x) <- paste0(library,"-", colnames(x))
  rownames(x) <- make.unique(rowData(x)$name)
  return(assays(x)[["z"]])
}

libs <- c("CD34_500_Day08", "CD34_500_Day14", "CD34_800_Day08", "CD34_800_Day14", "CD34_800_Day20")
SE <- do.call("cbind", lapply(libs, import_SE))
devs <- do.call("cbind", lapply(libs, import_CV)) %>% t()

stopifnot(all(rownames(devs) == colnames(SE)))

#-------------------------------------
# Everything below is code adapted from
# 
#-------------------------------------

#Binarize Sparse Matrix
binarizeMat <- function(mat){
  mat@x[mat@x > 0] <- 1
  mat
}

#LSI Adapted from fly-atac with information for re-projection analyses
calcLSI <- function(mat, nComponents = 50, binarize = TRUE, nFeatures = NULL, log = FALSE){
  
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

alignTrajectory <- function(df, trajectory, filter = 0.05, dof = 250, spar = 1){
  findClosest <- function(x, y, fitx, fity){
    distxy <- sqrt(rowSums(cbind((fitx - x)^2 + (fity - y)^2)))
    idxMin <- which.min(distxy)
    if(idxMin==1){
      idxMin <- idxMin + 1
    }else if(idxMin==length(fitx)){
      idxMin <- idxMin - 1
    }
    if(distxy[idxMin + 1]  < distxy[idxMin - 1]){
      diff <- 1
    }else{
      diff <- -1
    }
    data.frame(idx = idxMin, dist = distxy[idxMin], diff = diff)
  }
  dfAll <- data.frame()
  for(x in seq_along(trajectory)){
    #Subset
    dfx <- df[df$Group==trajectory[x],]
    #Mean Diff Filter
    xmean <- colMeans(dfx[,c(1,2)])
    diffx <- sqrt(colSums((t(dfx[,1:2]) - xmean)^2))
    dfx <- dfx[which(diffx <= quantile(diffx,1 - filter)),]
    #Get diff
    if(x!=length(trajectory)){
      xmean1 <- colMeans(df[df$Group==trajectory[x+1],c(1,2)])
      diffx1 <- sqrt(colSums((t(dfx[,1:2]) - xmean1)^2))
      dfx$time <- (1 - getQuantiles(diffx1)) + x
    }else{
      xmean1 <- colMeans(df[df$Group==trajectory[x-1],c(1,2)])
      diffx1 <- sqrt(colSums((t(dfx[,1:2]) - xmean1)^2))
      dfx$time <- getQuantiles(diffx1) + x
    }
    dfAll <- rbind(dfAll , dfx)
  }
  sx <- smooth.spline(dfAll$time, dfAll$x, df = dof, spar = spar)
  sy <- smooth.spline(dfAll$time, dfAll$y, df = dof, spar = spar)
  dfFit <- data.frame(x = sx[[2]], y = sy[[2]], t = seq_along(sy[[2]]))
  dfTrajectory <- df[df$Group %in% trajectory,]
  dfTime <- lapply(seq_len(nrow(dfTrajectory)), function(x){
    findClosest(dfTrajectory[x,1],dfTrajectory[x,2], dfFit[,1],dfFit[,2])
  }) %>% Reduce("rbind",.)
  dfTime$distQ <- getQuantiles(dfTime$dist)
  dfTrajectory$pseudotime <- 100*getQuantiles(dfTime$idx + matrixStats::rowProds(as.matrix(dfTime[,c("diff","distQ")])))
  
  out <- list(trajectory=dfTrajectory, fitTrajectory = dfFit)
}


projectLSI <- function(mat, lsi, log = FALSE){   
  
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

getQuantiles <- function(x){
  trunc(rank(x))/length(x)
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
nPCs1 <- 1:25
nPCs2 <- 1:25
resolution <- 0.8 #clustering resolution
nTop <- 25000 #number of variable peaks

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
clust2 <- louvainIgraphClusters(lsi2[[1]], 30)
length(unique(clust2))

#UMAP
set.seed(1)
umap <- umap::umap(
  lsi2$matSVD[,1:25], 
  n_neighbors = 55, 
  min_dist = 0.45, 
  metric = "euclidean", 
  verbose = TRUE    )
set.seed(10)

# Only multiplying by -1 from the umap output to be more consistent with a previous version
plot_df <- data.frame(umap$layout*-1, clusters = clust2, colData(SE),
                      data.matrix(devs[,c("GATA1", "CEBPA", "SPI1", "RUNX1", "PAX5", "KLF1", "IRF8")]))

# Plot TFs
p1 <- ggplot(plot_df, aes(x= X1, y = X2, color = ifelse(GATA1 > 10, 10, ifelse(GATA1 < -10, -10, GATA1)))) +
  geom_point_rast(raster.dpi = 500, size = 0.1) +
  scale_color_gradientn(colors = jdb_palette("brewer_spectra")) +
  theme_void() +theme(legend.position = "none") 

p2 <- ggplot(plot_df, aes(x= X1, y = X2, color = ifelse(SPI1 > 10, 10, ifelse(SPI1 < -10, -10, SPI1)))) +
  geom_point_rast(raster.dpi = 500, size = 0.1) +
  scale_color_gradientn(colors = jdb_palette("brewer_spectra")) +
  theme_void() +theme(legend.position = "none") 

p3 <- ggplot(plot_df, aes(x= X1, y = X2, color = ifelse(RUNX1 > 10, 10, ifelse(RUNX1 < -10, -10, RUNX1)))) +
  geom_point_rast(raster.dpi = 500, size = 0.1) +
  scale_color_gradientn(colors = jdb_palette("brewer_spectra")) +
  theme_void() +theme(legend.position = "none") 

p4 <- ggplot(plot_df, aes(x= X1, y = X2, color = ifelse(CEBPA > 10, 10, ifelse(CEBPA < -10, -10, CEBPA)))) +
  geom_point_rast(raster.dpi = 500, size = 0.1) +
  scale_color_gradientn(colors = jdb_palette("brewer_spectra")) +
  theme_void() +theme(legend.position = "none") 

cowplot::ggsave2(cowplot::plot_grid(p1, p2, p3, p4, nrow = 2), filename = "../plots/4TFs.png",
                width = 2.5, height = 2.5, units = "in", dpi = 1000)

re_anno <- c(
  'ery2', 'prog_my', 'ery4', 'ery6', 'my1',
  'my2','ery1','prog', 'prog_ery', 'ery3',
  'mix','my3','my2','my3','ery5',
  'my4'
)
names(re_anno) <- paste0("mc", as.character(1:16))
plot_df$new_name <- re_anno[as.character(plot_df$clusters)]

mean_df <- plot_df %>% group_by(clusters) %>%
  summarize(X1 = mean(X1), X2 = mean(X2))
ggplot(mean_df, aes(x= X1, y = X2, color = clusters, label = clusters)) +
  geom_text()
ggplot(plot_df, aes(x= X1, y = X2, color = clusters, label = clusters)) +
  geom_point() + geom_text(data = mean_df, color= "black")

mean_df2 <- plot_df %>% group_by(new_name) %>%
  summarize(X1 = mean(X1), X2 = mean(X2))
ggplot(mean_df2, aes(x= X1, y = X2, color = new_name, label = new_name)) +
  geom_text()

# Semi-supervised trajectory inference
trajectory1 <- c('prog', 'prog_ery', 'ery1', 'ery2', 'ery3', 'ery4', 'ery5', 'ery6')
trajectory2 <- c('prog', 'prog_my', 'my1', 'my2', 'my3', 'my4')

#Align single cells to Trajectory
df <- data.frame(row.names = rownames(plot_df),
                 x = plot_df$X1,
                 y = plot_df$X2,
                 Group = plot_df$new_name)
trajAligned1 <- alignTrajectory(df, trajectory1)
trajAligned2 <- alignTrajectory(df, trajectory2)

# simple function to add NAs to cells not in trajecotry
augment <- function(df_all, df_t){
  need <- !(rownames(df_all) %in% rownames(df_t))
  dfa <- df_all[need,c("X1", "X2", "new_name")]; dfa$ps <- NA
  colnames(dfa) <- c("x", "y", "Group", "pseudotime")
  total_df <- rbind(dfa, df_t)
  return(total_df)
}

df_ery <- augment(plot_df, trajAligned1[[1]])
df_my <- augment(plot_df, trajAligned2[[1]])

# Extract the pseudotime
dfT_ery <- trajAligned1[[2]]
dfT_my <- trajAligned2[[2]]

# Make plots for QC and for pseudotime
pB_name <- ggplot(plot_df, aes(x= X1, y = X2, color = clusters)) +
  geom_point_rast(raster.dpi = 500, size = 0.1) +
  theme_void() +theme(legend.position = "none") 

pB_mtdna <- ggplot(shuf(plot_df), aes(x= X1, y = X2, color = log10(mtDNAcoverage))) +
  geom_point_rast(raster.dpi = 500, size = 0.1) +
  scale_color_gradientn(colors = jdb_palette("solar_extra"), limits = c(0,3)) +
  theme_void() + theme(legend.position = "none") 

p_PS1 <- ggplot(df_ery, aes(x,y,color=pseudotime)) + 
  geom_point_rast(raster.dpi = 500, size = 0.1) +
  theme_void() +
  viridis::scale_color_viridis(na.value = "lightgrey") +
  geom_path(data=data.frame(dfT_ery), aes(x,y,color=NULL), size= 0.2, 
            arrow = arrow(type = "open", angle = 30, length = unit(0.05, "inches"))) +
  theme(legend.position = "none") 

p_PS2 <- ggplot(df_my, aes(x,y,color=pseudotime)) + 
  geom_point_rast(raster.dpi = 500, size = 0.1) +
  theme_void() +
  viridis::scale_color_viridis(na.value = "lightgrey") +
  geom_path(data=data.frame(dfT_my), aes(x,y,color=NULL), size= 0.2, 
            arrow = arrow(type = "open", angle = 30, length = unit(0.05, "inches"))) +
  theme(legend.position = "none")

cowplot::ggsave2(cowplot::plot_grid(pB_name, pB_mtdna, p_PS1, p_PS2, nrow = 2),
                filename = "../plots/trajectories_QC.png",
                width = 2.5, height = 2.5, units = "in", dpi = 1000)

# Export pseudotime inference for later use
save(plot_df, df_ery, df_my, dfT_ery, dfT_my,
     file = "../output/trajectory_inferences.18march2020.rda")

# Make final supplemental panels
plot_df$density <- get_density(plot_df$X1, plot_df$X2)
pdens <- ggplot(plot_df, aes(x= X1, y = X2, color = density)) +
  geom_point_rast(raster.dpi = 500, size = 0.1) +
  viridis::scale_color_viridis(na.value = "lightgrey") +
  theme_void() + theme(legend.position = "none") 

plib <- ggplot(shuf(plot_df), aes(x= X1, y = X2, color = library, label = library)) +
  geom_point_rast(raster.dpi = 500, size = 0.1) + theme_void() + 
  scale_color_manual(values = c("orange", "red","dodgerblue", "#008080", "purple4")) +
  theme(legend.position = "none")

cowplot::ggsave2(cowplot::plot_grid(pdens, plib, nrow = 1), filename = "../plots/twoSupp.png",
                width = 3.4, height = 1.7, units = "in", dpi = 1000)






