library(Signac)
library(Seurat)
library(data.table)
library(dplyr)
library(Matrix)
library(SummarizedExperiment)

# Function to import scATAC gene scores for label transfer
import_genescores_scATAC <- function(lib){
  mat <- readRDS(paste0("../../../mtscATACpaper_large_data_files/intermediate/",lib,".gene_activities.rds"))
  barcode_qc <- fread(paste0("../output/barcode_qc/",lib,".barcode_qc.tsv"))
  bcs <- barcode_qc %>% filter(keep) %>% pull(sample)
  mat_filt <- mat[,bcs]
  colnames(mat_filt) <- paste0(lib, "_", colnames(mat_filt))
  mat_filt
}

import_peaks_scATAC <- function(lib){
  mat <- readRDS(paste0("../../../mtscATACpaper_large_data_files/intermediate/",lib,".rds"))
  barcode_qc <- fread(paste0("../output/barcode_qc/",lib,".barcode_qc.tsv"))
  bcs <- barcode_qc %>% filter(keep) %>% pull(sample)
  mat_filt <- mat[,bcs]
  colnames(mat_filt) <- paste0(lib, "_", colnames(mat_filt))
  mat_filt
}

# Import RNA
pbmc.rna <- readRDS("../../../mtscATACpaper_large_data_files/intermediate/13March2020_recluster_reannotated_10xv3.rds")

# Import scATAC data
peaksmat <- SummarizedExperiment::cbind(import_peaks_scATAC("PBMC_H10"), import_peaks_scATAC("PBMC_H9"))
counts_mat <- assays(peaksmat)[["counts"]]
rownames(counts_mat) <- paste0("peak", as.character(1:dim(peaksmat)[1]))

# Prepare gene scores
h10 <- import_genescores_scATAC("PBMC_H10")
h9 <- import_genescores_scATAC("PBMC_H9")
genes_common <- intersect(rownames(h10), rownames(h9))

# Make Seurat object
pbmc.atac <- CreateSeuratObject(counts = counts_mat, assay = "ATAC", project = "10x_ATAC")
pbmc.atac[['ACTIVITY']] <- CreateAssayObject(counts = cbind(h9[genes_common,],h10[genes_common,]))
DefaultAssay(pbmc.atac) <- 'ATAC'

# Dimension reduction with signac
pbmc.atac <- RunTFIDF(pbmc.atac)
pbmc.atac <- FindTopFeatures(pbmc.atac, min.cutoff = 'q75')
pbmc.atac <- RunSVD(
  object = pbmc.atac,
  assay = 'ATAC',
  reduction.key = 'LSI_',
  reduction.name = 'lsi'
)
pbmc.atac <- RunUMAP(object = pbmc.atac, reduction = 'lsi', dims = 2:30)


# Do label transfer
transfer.anchors <- FindTransferAnchors(reference = pbmc.rna, query = pbmc.atac, features = VariableFeatures(object = pbmc.rna), 
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = pbmc.rna$celltype,
                                     weight.reduction = pbmc.atac[['lsi']],
                                     dims = 2:30)
pbmc.atac <- AddMetaData(pbmc.atac, metadata = celltype.predictions)
hist(pbmc.atac$prediction.score.max)
saveRDS(pbmc.atac@meta.data, file = "../output/PBMCatac_SignacSeurat_labelTransfer.rds")

DimPlot(
  object = pbmc.atac,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

FeaturePlot(
  object = pbmc.atac,
  features = "prediction.score.max"
)
