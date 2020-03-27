library(SummarizedExperiment)
library(Matrix)
library(dplyr)
library(BuenColors)
library(ComplexHeatmap)
library(circlize)
library(data.table)
library(ggrepel)
library(Seurat)
library(stringr)

source('../../global_functions/variant_calling.R')

"%ni%" <- Negate("%in%")
perc.rank <- function(x) trunc(rank(x))/length(x)

# Function to make cluster assignments from mtDNA heteroplasmy matrix
# Via seurat / cosine distance clustering
seuratSNN <- function(mat, resolution = 1, k.param = 10){ 
  set.seed(1)
  rownames(mat) <- make.unique(rownames(mat))
  obj <- FindNeighbors(mat, k.param = k.param, annoy.metric = "cosine")
  clusters <- FindClusters(object = obj$snn, resolution = resolution)
  return(as.character(clusters[,1]))
}

# Master function
process_patient_bcell_clone <- function(patient, resolution, knn){
  
  # Import and filter Bcell data
  mgatk_file <- paste0("../../../mtscATACpaper_large_data_files/source/mgatk_output/CLL_",patient,"_CD19pos_v12-mtMask_mgatk.rds")
  SE <- readRDS(mgatk_file)
  SE <- SE[,colData(SE)$depth >= 20]
  dim(SE)
  
  # Call mutations and determine abundance
  mut_se <- call_mutations_mgatk(SE, stabilize_variance = FALSE)
  misc_df <- data.frame(rowData(mut_se))
  filter_df <- misc_df %>%  filter(n_cells_conf_detected >= 5 & strand_correlation >= 0.65 & log10(vmr) > -2)
  dim(filter_df)
  
  p1 <- ggplot(misc_df %>% filter(n_cells_conf_detected >= 5), aes(x = strand_correlation, y = log10(vmr), color = log10(vmr) > -2 & strand_correlation > 0.65)) +
    geom_point(size = 0.4) + scale_color_manual(values = c("black", "firebrick")) +
    labs(color = "HQ", x = "Strand concordance", y = "log VMR") +
    pretty_plot(fontsize = 7) + L_border() +
    theme(legend.position = "bottom") + 
    geom_vline(xintercept = 0.6, linetype = 2) +
    geom_hline(yintercept = -2, linetype = 2) + theme(legend.position = "none")
  cowplot::ggsave(p1, file = paste0("../plots/CLL_",patient,"_CD19pos_mutation_viz.pdf"), width = 1.2, height = 1.2)
  
  afp <- as.matrix(data.matrix(assays(mut_se)[["allele_frequency"]][as.character(filter_df$variant),]))
  
  # Call clusters using the function above
  clusters <- seuratSNN(t(sqrt(afp)), resolution = resolution, k.param = knn) 
  table(clusters)
  clusters <- str_pad(as.character(as.numeric(factor(clusters))), 2, pad = "0")
  
  
  if(patient == "PT1"){
    var_order <- c("3179G>A",
                   "12213G>A", "4716C>A", "16519T>C",
                   "12067C>T","5809G>A",
                   "2702G>A", "11711G>A", 
                   "14858G>A","1260A>G","16390G>A",
                   "1872T>C","5140G>A","2587G>A","9950T>C","16052C>T","5835G>A","5274G>A"
                   
    )
  }
  
  if(patient == "PT2"){
    var_order <- c(
      "4853G>A", "12980G>A",  # one offs
      "11711G>A", "10367C>T","3496G>A", 
      "7775G>A",  "4429G>A", "4972G>A", "2702G>A","8020G>A", "14641A>G","12236G>A",
      "11982T>C", "1031G>A","15437G>A","3244G>A",   "2989G>A","1782G>A","16390G>A","5293G>A",
      "8858G>A" ,"13513G>A","16391G>A","2275T>C"
    )
    
    # Reorder clusters for presentation
    clusters <- clusterss
    nn <- c("a", "c", "g", "f", "b", "e", "d"); names(nn) <- sort(unique(clusters))
    clustersn <- nn[clusters]
    clusters <- clustersn
  }
  
  setdiff(var_order, rownames(afp))
  setdiff(rownames(afp), var_order)
  
  # Create pretty colors
  set.seed(2)
  names_clusters <- unique(sort((clusters)))
  vec_go <- c("grey", jdb_palettes[["flame_light"]], "dodgerblue2")[1:length(unique(clusters))]
  names(vec_go) <- names_clusters
  
  # Make data.frame for stuff
  df <- data.frame(
    cell_id = colnames(afp), 
    cluster_id = clusters
  ) %>% arrange(desc(clusters))
  
  ha_col <- HeatmapAnnotation(cell = as.character(df$cluster_id),
                              col = list(cell = vec_go))
  
  
  
  # Make heatmap
  afp[afp>0.1] <- 0.1
  afp[afp<0.01] <- 0.01
  
  png(paste0("../plots/",patient,"_scATAC_bcell_mtDNA_heatmap.png"), width=4, height=2, unit = "in", res = 500)
  hm <- Heatmap((data.matrix(afp[var_order,as.character(df$cell_id)])), #
                col=as.character(jdb_palette("solar_rojos",type="continuous")),
                show_row_names = TRUE, 
                top_annotation=ha_col,
                name = "AF",
                row_names_gp = gpar(fontsize = 4),
                cluster_rows = FALSE, 
                cluster_columns = FALSE,
                use_raster = FALSE,
                show_column_names = FALSE)
  hm 
  dev.off()
  
  pdf(paste0("../plots/",patient,"_scATAC_bcell_mtDNA_heatmap.pdf"), width=4, height=2)
  hm 
  dev.off()
  
  # Export data
  saveRDS(df, file = paste0("../output/",patient,"_clone_definition.rds"))
  
  mut_df_o <- data.frame(
    cluster = clusters, 
    t(data.matrix(assays(mut_se)[["allele_frequency"]][as.character(filter_df$variant),]))
  ) 
  saveRDS(mut_df_o, file = paste0("../output/",patient,"_specificVariants_forSupplement.rds"))
  
}

# Run function for two patients -- hyper parameter tweaks were determined by finding a 'reasonable' number of sub clones
# That seemed real by eye
process_patient_bcell_clone("PT1", 0.2, 20) # Default gives thousands of 'sub clones'
process_patient_bcell_clone("PT2", 1.0, 30) # Default gives ~3 subclones


