library(SummarizedExperiment)
library(chromVAR)
library(BSgenome.Hsapiens.UCSC.hg19)
library(chromVARmotifs)
library(motifmatchr)
source("/data/aryee/caleb/tenx-scatac/paper/mtscATACpaper_reproducibility/global_functions/helper_SE_functions.R")
data("human_pwms_v2")

# Basic function to import summarized experiment
import_SE <- function(library){
  x <- readRDS(paste0("../../../mtscATACpaper_large_data_files/intermediate/",library,".rds"))
  colData(x)$library <- library
  colnames(x) <- paste0(library,"-", colnames(x))
  return(x)
}


# Import all of the in vitro heme data
SE <- do.call("cbind", lapply(c("CD34_500_Day08", "CD34_500_Day14", "CD34_800_Day08", "CD34_800_Day14", "CD34_800_Day20"), import_SE))

# Start doing chromVAR things
SE <- addGCBias(SE, BSgenome.Hsapiens.UCSC.hg19)
counts <- SummarizedExperiment::assays(SE)[["counts"]]

# Produce matrices and export
bg_peaks <- background_peaks_CL(Matrix::rowSums(counts) + 1, mcols(rowRanges(SE))$bias, niterations = 50)
saveRDS(bg_peaks, file = "../output/chromVAR/chromvar_backgroundpeaks_invitro.rds")

mm <- matchMotifs(human_pwms_v2, SE, BSgenome.Hsapiens.UCSC.hg19)
saveRDS(mm, file = "../output/chromVAR/motifmatches_humanPWMS2.rds")

