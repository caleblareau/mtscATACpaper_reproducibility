library(chromVAR)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)
library(tools)
library(data.table)
library(dplyr)
register(MulticoreParam(4)) # Use 4 cores

args <- commandArgs(trailingOnly = TRUE)
i <- as.numeric(args[1])

"%ni%" <- Negate("%in%")
file <- list.files("rds_SE", full.names = TRUE)[i]
SE <- readRDS(file)

bg <- readRDS("chromVAR/backgroundpeaks_ery100.rds")

# Do it for a specific experiment
outname <- gsub(".rds", "", list.files("rds_SE", full.names = FALSE)[i])

# Gently modify counts to work with chromVAR
s <- Matrix::summary(counts(SE))
new <- which(1:dim(SE)[1] %ni% s$i)
SE_new <- SummarizedExperiment(
  rowRanges = rowRanges(SE),
  assays = list(counts =   Matrix::sparseMatrix(
    i = c(s$i, new, 451283),
    j = c(s$j, sample(1:dim(SE)[2], length(new), replace = TRUE), 2), 
    x = c(ifelse(s$x > 3, 3, s$x), rep(1, length(new)), 1)     
  )),
  colData = SE@colData
)
ss <- Matrix::summary(counts(SE_new))
new <- which(1:dim(SE_new)[1] %ni% ss$i)

mm <- readRDS("chromVAR/motifmatches_humanPWMS2.rds")
dev <- computeDeviations(SE_new, mm, background_peaks = bg, expectation = computeExpectations(SE_new))
saveRDS(dev, file = paste0("chromVAR/", outname, "_tf_deviations.rds"))