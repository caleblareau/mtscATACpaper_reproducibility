library(Matrix)
library(SummarizedExperiment)
library(tidyverse)
library(umap)
library(edgeR)
library(FNN)
library(matrixStats)
library(igraph)
set.seed(1)
load("15November2019_projection_5k_CLL.rda")

#annovec <- c("CD4memory", "CD4activated", "pDC", "CD4naive", "DC", "Bcell", "CD16mono", "CD14mono", "NKcell", "CD8tcell")
#names(annovec) <- paste0("mc", as.character(1:10))


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


classify_from_reference <- function(A, B){
  
  # Get the pair-wise euclidean distance msot-alike
  euklDist <- t(sqrt(apply(array(apply(B,1,function(x){(x-t(A))^2}),c(ncol(A),nrow(A),nrow(B))),2:3,sum)))
  colnames(euklDist) <- rownames(A)
  colnames(euklDist)[max.col(-1*euklDist, 'first')]  -> vec
  return(vec)
}

import_classify <- function(se_file, outfile){
  seDisease <- readRDS(se_file)
  matProjectLSI <- assay(seDisease[lsiPeaks,])
  
  #LSI Project + classify
  lsiProjection <- projectLSI(matProjectLSI, lsiReference)
  projected_clustersbasic <- classify_from_reference(t(means_lsi)[,1:25],data.matrix(lsiProjection)[,1:25])
  projected_clusters <- annovec[as.character(projected_clustersbasic)]
  
  #UMAP Projection
  #Set Seed Prior to umap_transform (see uwot github)
  set.seed(1)
  umapProjection <- round(predict(umap, data.matrix(lsiProjection[,1:25])), 2)
  odf <- data.frame(barcode = rownames(umapProjection), umapProjection, cluster = projected_clusters)
  names(odf) <- c("barcode", "umap1", "umap2", "cluster")

  write.table(odf, paste0("projection_dfs/", outfile, "_final.tsv"), 
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  outfile
}
import_classify("rds_SE/Patient1_CLL_15min_QC.rds", "Pt1_15min")
import_classify("rds_SE/Patient1_CLL_5min_QC.rds", "Pt1_5min")
import_classify("rds_SE/Patient1_CLL_neg_QC.rds", "Pt1_CD19neg")
import_classify("rds_SE/Patient1_CLL_pos_QC.rds", "Pt1_CD19pos")

import_classify("rds_SE/Patient2_CLL_neg_QC.rds", "Pt2_CD19neg")
import_classify("rds_SE/Patient2_CLL_pos_QC.rds", "Pt2_CD19pos")


