library(Signac)
library(Seurat)
library(ComplexHeatmap)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(BuenColors) # For caleb Plotting

"%ni%" <- Negate("%in%")
source("../../global_functions/variant_calling.R") 

# Function to import the peak x cells count matrix that's been gently modified (gzip) from 
# CellRanger-ATAC
import_atac_counts_gz <- function(directory){
  
  # Import cells and annotate counts matrix
  barcodes <- fread(paste0(directory,"/barcodes.tsv.gz"), header = FALSE)[[1]]
  mtx <- fread(paste0(directory,"/matrix.mtx.gz"), header = FALSE, skip = 3)
  counts <- sparseMatrix(i = mtx[[1]], j = mtx[[2]], x = mtx[[3]])
  colnames(counts) <- barcodes
  
  # Import peaks and annotate matrix
  drr <- fread(paste0(directory,"/peaks.bed.gz"), header = FALSE)
  rownames(counts) <- paste0(drr[["V1"]], ":", drr[["V2"]], "-", drr[["V3"]])
  counts
}

# Import everything
counts <- import_atac_counts_gz("../data/CRC_filtered_peak_bc_matrix/")
metadata <- read.csv(
  file = "../data/singlecell_sumstats/CRC_singlecell.csv.gz",
  header = TRUE,
  row.names = 1
)

# Make a seurat object with the inputs
crc <- CreateSeuratObject(
  counts = counts,
  assay = 'peaks',
  project = 'ATAC',
  min.cells = 1,
  meta.data = metadata
)

# Augment QC
crc$pct_reads_in_peaks <- crc$peak_region_fragments / crc$passed_filters * 100
crc$pct_reads_in_DNase <- crc$DNase_sensitive_region_fragments / crc$passed_filters * 100

# There will be a population of high fragment, low frip cells that we will choose to keep to examine downstream
qplot(log10(crc@meta.data$passed_filters), crc@meta.data$pct_reads_in_DNase) +
  pretty_plot() + L_border() + labs(x = "# Unique Fragments", y = "% Reads in DNase Hypersensitivity Sites") +
  geom_vline(xintercept = 3, color = "firebrick") +
  geom_hline(yintercept = 40, color = "firebrick")

crc_filt <- subset(crc, subset = passed_filters >= 1000 & pct_reads_in_DNase >= 40)
crc_filt

# Now subset based on mtDNA coverage for cells
mgatk_se <- readRDS("../../../mtscATACpaper_large_data_files/source/mgatk_output/CRC_v12-mtMask_mgatk.rds")

# Subset to HQ cells that exist so far
mgatk_se <- mgatk_se[,colnames(crc_filt)]

# Threshold based on abundance of mtDNA
crc_filt$mtDNA_depth <- mgatk_se$depth

# Could probably replace this with a violin plot
qplot(crc_filt$mtDNA_depth) + scale_x_log10(breaks = c(1,10,100)) +
  geom_vline(xintercept = 10, color = "firebrick") +
  pretty_plot() + L_border() + labs(x = "mean mtDNA depth/bp", y = "count")

# Final filter -- mgatk, Signac, and fragments all
crc_filt2 <- subset(crc_filt, mtDNA_depth >=10)
write.table(data.frame(colnames(crc_filt2)), 
            file = "../output/CRC_filtered_barcodes.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Rename for convenience
crc <- crc_filt2

# Run semi-standard Signac dimension reduction, clustering, etc. 
crc <- RunTFIDF(crc)
crc <- FindTopFeatures(crc, min.cutoff = 'q50')
crc <- RunSVD(
  object = crc,
  assay = 'peaks',
  reduction.key = 'LSI_',
  reduction.name = 'lsi',
  irlba.work=500
)

crc <- RunUMAP(object = crc, reduction = 'lsi', dims = 2:15)
crc <- FindNeighbors(object = crc, reduction = 'lsi', dims = 2:15)
crc <- FindClusters(object = crc, verbose = FALSE, resolution = 0.3)
DimPlot(object = crc, label = TRUE) + NoLegend() +
  scale_color_manual(values = c("dodgerblue3",  "forestgreen", "orange3", "firebrick", "purple3", "magenta"))

crc$log10_mtDNA_depth <- log10(crc$mtDNA_depth)
FeaturePlot(crc, features = c("pct_reads_in_DNase"))
FeaturePlot(crc, features = c("log10_mtDNA_depth"))

#-----------------------------
# Process gene activity scores
#-----------------------------

ga_rds <- "../../../mtscATACpaper_large_data_files/intermediate/CRC-geneActivites.rds"
if(!file.exists(ga_rds)){
  # extract gene coordinates from Ensembl, and ensure name formatting is consistent with Seurat object 
  gene.coords <- genes(EnsDb.Hsapiens.v75, filter = ~ gene_biotype == "protein_coding")
  seqlevelsStyle(gene.coords) <- 'UCSC'
  genebody.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')
  genebodyandpromoter.coords <- Extend(x = gene.coords, upstream = 2000, downstream = 0)
  
  # create a gene by cell matrix
  gene.activities <- FeatureMatrix(
    fragments = "../../../mtscATACpaper_large_data_files/source/cellranger_output/CRC_v12-mtMask_fragments.tsv.gz",
    features = genebodyandpromoter.coords,
    cells = colnames(crc),
    chunk = 20
  )
  
  # convert rownames from chromsomal coordinates into gene names
  gene.key <- genebodyandpromoter.coords$gene_name
  names(gene.key) <- GRangesToString(grange = genebodyandpromoter.coords)
  rownames(gene.activities) <- gene.key[rownames(gene.activities)]
  saveRDS(gene.activities, ga_rds)
} else{
  gene.activities <- readRDS(ga_rds)
}

# Put into Seurat object
crc[['RNA']] <- CreateAssayObject(counts = gene.activities)
crc <- NormalizeData(
  object = crc,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(crc$nCount_RNA)
)

DefaultAssay(crc) <- 'RNA'

FeaturePlot(
  object = crc,
  features = c("IL3RA", "HDC", "KIAA0226L", "KRT8", "RORA", "IL21R"),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)


#--- 
# Now call mtDNA variants
#--- 
mut_se <- call_mutations_mgatk(mgatk_se[,colnames(crc)])
misc_df <- data.frame(rowData(mut_se))
filter_df <- misc_df %>% dplyr::filter(n_cells_detected >= 5 & strand_correlation >= 0.65 & log10(vmr) > -2)
dim(filter_df)

full_df <- data.frame(
  crc@meta.data,
  crc@reductions$umap@cell.embeddings,
  t(assays(mut_se)[["allele_frequency"]][as.character(filter_df$variant),]),
  t(crc@assays$RNA@data[c('TREM1', 'EPCAM', "PTPRC", "IL1RL1", "KIT", "GATA3", "IL3RA", "HDC", "KIAA0226L", "KRT8", "RORA", "IL21R"),])
)

saveRDS(full_df, file = "../output/21March2020_signac_process.rds")
