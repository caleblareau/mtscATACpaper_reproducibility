options(warn=-1)
library(FNN)
library(umap)
library(SummarizedExperiment)
library(irlba)
library(BuenColors)
library(chromVAR)
library(BSgenome.Hsapiens.UCSC.hg19)

suppressMessages(suppressWarnings(library(tools)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(dplyr)))


# Function to import chromVAR (CV) Z scores into a more simple object for cbinding
importCV <- function(file){
  
  # Import and filter based on obvious peripherals
  CVdev <- readRDS(file)
  
  SE <- SummarizedExperiment(
    assays = list("z" = assays(CVdev)[["z"]]),
    rowData = rowData(CVdev),
    colData = colData(CVdev)
  )
  return(SE)
  
}


substrRight <- function(x, n) substr(x, nchar(x)-n+1, nchar(x))

# Function that creates background peak indicies via chromVAR
# But more convenient
background_peaks_CL<- function(fragments_per_peak, 
                               bias, 
                               niterations = 250, 
                               w = 0.1, 
                               bs = 50) {
  
  intensity <- log10(fragments_per_peak + 1)
  norm_mat <- matrix(c(intensity, bias), ncol = 2, byrow = FALSE)
  
  chol_cov_mat <- chol(cov(norm_mat))
  trans_norm_mat <- t(forwardsolve(t(chol_cov_mat), t(norm_mat)))
  
  # make bins
  bins1 <- seq(min(trans_norm_mat[, 1]), max(trans_norm_mat[, 1]), 
               length.out = bs)
  bins2 <- seq(min(trans_norm_mat[, 2]), max(trans_norm_mat[, 2]), 
               length.out = bs)
  
  bin_data <- do.call(rbind, lapply(seq_len(bs), 
                                    function(x) matrix(c(rep(bins1[x], bs), 
                                                         bins2), ncol = 2, 
                                                       byrow = FALSE)))
  
  bin_dist <- chromVAR:::euc_dist(bin_data)
  bin_p <- dnorm(bin_dist, 0, w)
  
  bin_membership <- nabor::knn(bin_data, query = trans_norm_mat, k = 1)$nn.idx
  
  bin_density <- chromVAR:::tabulate2(bin_membership, min_val = 1, max_val = bs^2)
  
  background_peaks <- chromVAR:::bg_sample_helper(bin_membership - 1, bin_p, bin_density, 
                                                  niterations)
  
  return(background_peaks)
  
}
