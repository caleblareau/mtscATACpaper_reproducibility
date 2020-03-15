library(Seurat)
library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
library(BuenColors)

set.seed(1)
# Script gently adapted from https://github.com/satijalab/Integration2019/blob/master/preprocessing_scripts/pbmc_10k_v3.R

import_scRNAseq <- function(dir_base, name = "rna"){
  data.dir <- paste0(dir_base)
  raw <- Read10X(data.dir = data.dir); colnames(raw) <- paste0(name, "-", colnames(raw))
  
  # import scrublet results
  singlets <- fread(paste0(dir_base, "/scrublet_out.tsv")) %>%
    data.frame() %>% dplyr::filter(!called) %>% pull(barcode)
  singlets_rename <-  paste0(name, "-", substr(singlets, 1, 16))
  raw <- raw[,singlets_rename]
  
  # Filter for singlet genes adn non-mitos
  rna <- CreateSeuratObject(counts = raw,  min.cells = 3, min.features = 200)

  # Add some more QC
  mito.features <- grep(pattern = "^MT-", x = rownames(x = rna), value = TRUE)
  percent.mito <- Matrix::colSums(x = GetAssayData(object = rna, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = rna, slot = 'counts'))
  rna$percent.mito <- percent.mito
  rna <- subset(x = rna, subset = nCount_RNA > 2000 & nCount_RNA < 20000 & percent.mito < 0.2)
  rna
  
}

# Run the function
rna <- import_scRNAseq("pbmc_10k_v3")

# QC
rna <- NormalizeData(rna)
rna <- FindVariableFeatures(rna, nfeatures = 3000)
rna <- ScaleData(rna)
rna <- RunPCA(rna, npcs = 100)
rna <- RunUMAP(rna, dims = 1:30)
rna <- FindNeighbors(rna, dims = 1:30)
rna <- FindClusters(rna)

DimPlot(rna, label  = TRUE) + scale_color_manual(values = jdb_palette("corona"))


# Visualize relevant features

# CL feature plot function
fpcl <- function(feature){
  FeaturePlot(rna, features = c(feature), min.cutoff = "q8", pt.size = 0.05) + 
    theme_void() +theme(legend.position = "none") 
}
ggsave(plot_grid(fpcl("CD14"), fpcl("FCGR3A"), fpcl("CD1C"),  fpcl("IL3RA"), fpcl("MX1"), ncol = 5), 
       file = "plots/Monocyte.png",width = 14, height = 3, dpi = 1000)

ggsave(plot_grid(fpcl("PPBP"),  fpcl("MS4A1"), fpcl("CD27"), fpcl("TCL1A"), fpcl("CD34"), ncol = 5), 
       file = "plots/rare_Bcell.png", width = 14, height = 3, dpi = 500)

ggsave(plot_grid(fpcl("CCL5"),  fpcl("NCR3"), fpcl("GZMK"), fpcl("GNLY"), fpcl("CD8A"), ncol = 5), 
       file = "plots/NKcell.png", width = 14, height = 3, dpi = 500)

ggsave(plot_grid(fpcl("CCR7"),  fpcl("CD4"), fpcl("S100A4"), fpcl("TRGC2"), fpcl("FOXP3"), ncol = 5), 
       file = "plots/Tcell.png", width = 14, height = 3, dpi = 500)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
rna <- CellCycleScoring(rna, s.features = s.genes, g2m.features = g2m.genes)

FeaturePlot(rna, features = "G2M.Score")
ggsave(plot_grid(fpcl("G2M.Score"), ncol = 1), 
       file = "../plots/G2M_score.png", width = 2.8, height = 3, dpi = 500)


FindMarkers(rna, ident.1 = "16", min.pct = 0.3) %>% head(25)
FindMarkers(rna, ident.1 = "5", ident.2 = "6", min.pct = 0.3) %>% head(25)

new.cluster.ids <- c(
  '0'='Memory_CD4_Tcell',
  '1'='CD14_monocyte',
  '2'='Naive_CD4_Tcell',
  '3'='CD14_monocyte',
  '4'='activated_Bcell',
  '5'='Cytotoxic_CD8_Tcell',
  '6'='GammaDelta_Tcell',
  '7'='NK_cell',
  '8'='Memory_Bcell',
  '9'="Dendritic_cell",
  '10'='Naive_CD4_Tcell',
  '11'='CD16_monocyte',
  '12'='Treg',
  '13'='IFNactive_monocyte',
  '14'='Platelet',
  '15'='pDC',
  '16'='Progenitor'
)

table(rna@meta.data$seurat_clusters)

Idents(rna) <- new.cluster.ids[as.character(rna@meta.data$seurat_clusters)]
rna$celltype <- Idents(rna)
table(rna@meta.data$celltype)

DimPlot(rna, label  = TRUE) + scale_color_manual(values = jdb_palette("corona"))
saveRDS(rna, "13March2020_recluster_reannotated_10xv3.rds")
