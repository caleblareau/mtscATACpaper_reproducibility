library(BuenColors)
library(data.table)
library(circlize)
library(ComplexHeatmap)
library(dplyr)
library(SummarizedExperiment)
library(Matrix)

# Import data
sdf <- readRDS("../output/21March2020_signac_process.rds")
se <- readRDS("../output/27March_mutSE_filt_CRC.rds")
afp <- assays(se)[["allele_frequency"]]

# set colors
names_clusters <- sort(unique(as.character(sdf$seurat_clusters)))
colors <- c("dodgerblue3", "purple4", "forestgreen", "orange3", "firebrick", "black")
names(colors) <- sort(names_clusters)

# Make data.frame for visualization
df <- data.frame(
  cell_id = colnames(afp), 
  cluster_id = as.character(sdf$seurat_clusters)
) 

ha_col <- HeatmapAnnotation(cell = as.character(df$cluster_id),
                            col = list(cell = colors))

afp[afp < 0.01] <- 0
afp[afp > 0.1] <- 0.1

pdf(paste0("../plots/complete_mito_heatmap.pdf"), width=3.2, height=2)
hm <- Heatmap((data.matrix(afp)[,as.character(df$cell_id)]),  
              col=as.character(jdb_palette("solar_rojos",type="continuous")),
              show_row_names = TRUE, 
              top_annotation=ha_col,
              cluster_columns = TRUE,
              name = "AF",use_raster = FALSE,
              row_names_gp = gpar(fontsize = 4),
              cluster_rows = TRUE, 
              show_column_names = FALSE)
hm 
dev.off()

