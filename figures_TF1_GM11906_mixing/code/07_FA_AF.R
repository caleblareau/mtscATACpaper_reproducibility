library(SummarizedExperiment)
library(data.table)
library(Matrix)
library(dplyr)
source("02a_assign_pull_function.R")

# Compute bulk vector of heteroplasmy given a summarized experiment
compute_bulk_vector <- function(SE){
  # Determinie key coverage statistics every which way
  ref_allele <- toupper(as.character(rowRanges(SE)$refAllele))
  adf <- assign_pull(SE)
  
  # Filter for relatively pure cells as we are testing the effect of FA
  cids <- adf %>% filter(gm11906_mean <= 0.002 & assign == "TF1") %>% pull(cell_id)
  SE <- SE[,colnames(SE) %in% cids]
  cov <- assays(SE)[["coverage"]]
  
  pos <- start(rowRanges(SE))
  getBulk <- function(letter){
    boo <- ref_allele != letter & ref_allele != "N"
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


# Visualize to confirm FA treatment does not induce false positive mutations
cor(pldf[,c(1:3)])

p1 <- ggplot(pldf, aes(x = noFA, y = FA01)) +
  geom_point_rast(size = 0.1, raster.dpi = 500) +
  pretty_plot(fontsize = 8) + L_border() 

cowplot::ggsave2(p1, file = "../plots/compare_noFA_FA01.pdf", width = 2, height = 1.8)

p2 <- ggplot(pldf, aes(x = noFA, y = FA1)) +
  geom_point_rast(size = 0.1, raster.dpi = 500) +
  pretty_plot(fontsize = 8) + L_border() 

cowplot::ggsave2(p2, file = "../plots/compare_noFA_FA1.pdf", width = 2, height = 1.8)

