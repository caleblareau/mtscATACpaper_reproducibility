library(chromVAR)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)
library(tools)
library(data.table)
library(dplyr)
register(MulticoreParam(4)) # Use 4 cores

# Pull i from the environment -- executed via a cluster
args <- commandArgs(trailingOnly = TRUE)
i <- as.numeric(args[1])

"%ni%" <- Negate("%in%")
samples <- c("CD34_500_Day08", "CD34_500_Day14", "CD34_800_Day08", "CD34_800_Day14", "CD34_800_Day20")

library <- samples[i]
file <- paste0("../../../mtscATACpaper_large_data_files/intermediate/",library,".rds")
SE <- readRDS(file)

# Import background peaks
bg <- readRDS("../output/chromVAR/chromvar_backgroundpeaks_invitro.rds")

# Gently modify counts to work with chromVAR -- make sure that each peak has a read; make sure no one peak has a crazy # of counts
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

mm <- readRDS("../output/chromVAR/motifmatches_humanPWMS2.rds")
dev <- computeDeviations(SE_new, mm, background_peaks = bg, expectation = computeExpectations(SE_new))
saveRDS(dev, file = paste0("../output/chromVAR/", library, "_tf_deviations.rds"))

