library(BuenColors)
library(circlize)
library(ComplexHeatmap)
library(SummarizedExperiment)

#---
# Make heatmap of variants
#---
mut_se <- readRDS("../output/27March_mutSE_filt_CRC.rds")
sdf <- readRDS("../output/21March2020_signac_process.rds")

af = data.matrix(assays(mut_se)[["allele_frequency"]])[,sdf$seurat_clusters == 0]
af[af>0.1] <- 0.1
muts_plot <- c(
  "16147C>T", "12889G>A", "9728C>T", "1227G>A", "6081G>A", "824T>C"
)
afp <- af[muts_plot,]
wrs <- afp[1,]*600000 + afp[2,]*50000 + afp[3,]*4000  + afp[4,]*300 + afp[5,]*20 + afp[5,]*5
names(wrs) <- colnames(afp)
swrs <- sort(wrs, decreasing = TRUE)

pdf(paste0("../plots/heatmap_mtDNAvars.pdf"), width=3, height=1.5)
hm <- Heatmap(afp[,names(swrs)], 
              col=as.character(jdb_palette("solar_rojos",type="continuous")),
              show_row_names = TRUE, 
              cluster_columns = FALSE,
              cluster_rows = FALSE,
              row_names_gp = gpar(fontsize = 5),
              column_names_gp = gpar(fontsize = 5),
              show_column_names = FALSE)
hm
dev.off()


