library(BuenColors)
library(ggrastr)
library(dplyr)
library(data.table)

source("../../global_functions/get_allele_freq_mat.R")

import_munge <- function(exp){
  
  mgatk_file <- paste0("../../../mtscATACpaper_large_data_files/source/mgatk_output/", exp, "_v12-mtMask_mgatk.rds")
  proj_file <- paste0("../output/projection_dfs/", exp, "_projection.tsv")
  
  se <- readRDS(mgatk_file)
  af <- computeAFMutMatrix(se)
  
  bcsdf <- fread(proj_file, header = TRUE) %>% data.frame()
  
  # Annotate with mtDNA information
  df2 <- data.frame(barcode = colnames(af),
                    data.matrix(t(af[c("5140G>A", "1260A>G", "14858G>A", "1872T>C"),] >= 0.2)), depth = colData(se)$depth) %>% filter(depth > 20)
  mdfff <- merge(bcsdf, df2, by.x = "barcode", by.y = "barcode")
  return(mdfff)
}

cd19neg <- import_munge("CLL_PT1_bulk1")
bulk1 <- import_munge("CLL_PT1_bulk2")
bulk2 <- import_munge("CLL_PT1_CD19neg")

# Combine samples for visualization
df <- rbind(cd19neg, bulk1, bulk2)
df$score <- ifelse(df$X5140G.A, "5140G>A", ifelse(df$X1260A.G, "1260A>G", ifelse(df$X14858G.A, "14858G>A", ifelse(df$X1872T.C, "1872T>C", ""))))

p1a <- ggplot(rbind(df[df$score == "",], shuf(df[df$score != "",])), aes(x = umap1, y =umap2, color = score)) +
  geom_point_rast(size = 0.2, raster.dpi = 1000) + pretty_plot(fontsize = 7) + L_border() +
  scale_color_manual(values = c("lightgrey", "dodgerblue4", "dodgerblue", "orange2","firebrick")) + 
  theme(legend.position = "none") + labs(x = "UMAP1", y = "UMAP2")
cowplot::ggsave2(p1a, file = "../plots/UMAP_CD19neg_CLL1_projection.pdf", 
                width = 1.7, height = 1.7)

