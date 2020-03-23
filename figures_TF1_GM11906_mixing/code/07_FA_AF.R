library(SummarizedExperiment)
library(data.table)
library(Matrix)
library(dplyr)

# Compute bulk vector of heteroplasmy given a summarized experiment
compute_bulk_vector <- function(SE){
  # Determinie key coverage statistics every which way
  cov <- assays(SE)[["coverage"]]
  ref_allele <- toupper(as.character(rowRanges(SE)$refAllele))
  
  getBulk <- function(letter){
    boo <- ref_allele != letter & ref_allele != "N"
    pos <- start(rowRanges(SE))
    variant_name <- paste0(as.character(pos), ref_allele, ">", letter)[boo]
    vec <- (Matrix::rowSums(assays(SE)[[paste0(letter, "_counts_fw")]] + assays(SE)[[paste0(letter, "_counts_rev")]]) / Matrix::rowSums(cov))[boo]
    names(vec) <- variant_name
    return(vec)
  }
  c(getBulk("A"), getBulk("C"), getBulk("G"), getBulk("T"))
}

noFA <- compute_bulk_vector(readRDS("../../../mtscATACpaper_large_data_files/source/mgatk_output/Mix_1_CR-mtMask_mgatk.rds"))
FA01 <- compute_bulk_vector(readRDS("../../../mtscATACpaper_large_data_files/source/mgatk_output/Mix_3_Fix_12h_CR-mtMask_mgatk.rds"))
FA1 <- compute_bulk_vector(readRDS("../../../mtscATACpaper_large_data_files/source/mgatk_output/Mix_4_Fix_12h_CR-mtMask_mgatk.rds"))

pldf <- data.frame(
  noFA,
  FA01,
  FA1,
  variant = names(noFA)
)

cor(pldf[,c(1:3)])

ggplot(pldf, aes(x = noFA, y = FA01)) +
  geom_point()
