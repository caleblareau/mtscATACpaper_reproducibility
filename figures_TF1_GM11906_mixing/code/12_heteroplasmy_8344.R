library(SummarizedExperiment)
library(BuenColors)
library(dplyr)
library(Matrix)
library(ggrastr)

perc.rank <- function(x) trunc(rank(x))/length(x)

# Import mtscATAC summarized experiment data w/ meta data in the colData
SE <- readRDS("../../../mtscATACpaper_large_data_files/intermediate/GM11906_combinedSE.rds")

mscatac <- data.frame(
  af = colData(SE)$merrf_af,
  cov = colData(SE)$merrf_cov,
  assay = "mscATAC"
)

mscatac$coverage_quantile <- perc.rank(mscatac$cov)

# Functiont to pull the heteroplasmy from the C1
pull_het_SE <- function(SE, assay = "C1"){
  cov <- assays(SE)[["coverage"]]+ 0.001
  ref_allele <- as.character(rowRanges(SE)$refAllele)
  
  getMutMatrix <- function(letter){
    mat <- (assays(SE)[[paste0(letter, "_counts_fw")]] + assays(SE)[[paste0(letter, "_counts_rev")]]) / cov
    rownames(mat) <- paste0(as.character(1:dim(mat)[1]), toupper(ref_allele), ">", letter)
    return(mat[toupper(ref_allele) != letter,])
  }
  
  mmat <- rbind(getMutMatrix("A"), getMutMatrix("C"), getMutMatrix("G"), getMutMatrix("T"))
  
  data.frame(af = mmat["8344A>G",colMeans(cov) > 20], cov = round(cov[8344,colMeans(cov) > 20], 1), assay = assay)
}

C1het <- pull_het_SE(readRDS("../data/GM11906_C1.rds"), "C1")
C1het$coverage_quantile <- perc.rank(C1het$cov)

# Import FISH data
FISH <- read.csv("../data/heteroplasmy_and_coverage_fov09.csv", col.names = c("af", "cov"))
FISH$assay <- "FISH"
FISH$coverage_quantile <- perc.rank(FISH$cov)

# Determine medians
plot_df <- rbind(mscatac, C1het, FISH)
median_df <- plot_df %>% group_by(assay) %>%
  summarize(af = mean(af))

# Make a plot
p1 <- ggplot(plot_df, aes(x = assay, y = af*100, col = coverage_quantile)) +
  geom_quasirandom_rast(size = 0.2, raster.dpi = 1000) + pretty_plot(fontsize = 8) + L_border() + 
  geom_point(data = median_df, inherit.aes = FALSE, aes(x = assay, y = af*100), color = "black", shape = 124, size = 4) +
  labs(x = "", y = "m.8344A>G Heteroplasmy (%)") +
  theme(axis.text.y = element_text(angle = 45, hjust = 1)) + 
  theme(legend.position = "none") + coord_flip() + 
  scale_color_gradientn(colors = jdb_palette("solar_basic")) +
  geom_hline(yintercept = 44, linetype = 2)

cowplot::ggsave2(p1, file = "../plots/MERRFallele_quant_methods.pdf", width = 2.7, height = 1.5)
