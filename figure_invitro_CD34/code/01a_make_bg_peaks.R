library(SummarizedExperiment)
library(chromVAR)
library(BSgenome.Hsapiens.UCSC.hg19)
library(chromVARmotifs)
library(motifmatchr)

data("human_pwms_v2")
import_SE <- function(library){
  x <- readRDS(paste0("rds_SE/",library,".rds"))
  colData(x)$library <- library
  colnames(x) <- paste0(library,"-", colnames(x))
  return(x)
}

SE <- do.call("cbind", lapply(c("CD34_500", "CD34_500_Epo", "CD34_800", "CD34_800_CC100", "CD34_800_Epo"), import_SE))

SE <- addGCBias(SE, BSgenome.Hsapiens.UCSC.hg19)
counts <- SummarizedExperiment::assays(SE)[["counts"]]

source("/data/aryee/caleb/biorad/all_heme/helper_fns.R")
bg_peaks <- background_peaks_CL(Matrix::rowSums(counts) + 1, mcols(rowRanges(SE))$bias, niterations = 100)
saveRDS(bg_peaks, file = "chromVAR/backgroundpeaks_ery100.rds")

mm <- matchMotifs(human_pwms_v2, SE, BSgenome.Hsapiens.UCSC.hg19)
saveRDS(mm, file = "chromVAR/motifmatches_humanPWMS2.rds")

