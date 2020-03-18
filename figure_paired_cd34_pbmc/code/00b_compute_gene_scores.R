library(Signac)
library(Seurat)

library(data.table)
library(SummarizedExperiment)

if(FALSE){
  library(GenomeInfoDb)
  library(EnsDb.Hsapiens.v75)
  
  # extract gene coordinates from Ensembl, and ensure name formatting is consistent with Seurat object 
  gene.coords <- genes(EnsDb.Hsapiens.v75, filter = ~ gene_biotype == "protein_coding")
  seqlevelsStyle(gene.coords) <- 'UCSC'
  genebody.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')
  genebodyandpromoter.coords <- Extend(x = gene.coords, upstream = 2000, downstream = 0)
  saveRDS(genebodyandpromoter.coords, "../data/genebodyandpromoter.coords_EnsDb.Hsapiens.v75.rds")
}

genebodyandpromoter.coords <- readRDS("../data/genebodyandpromoter.coords_EnsDb.Hsapiens.v75.rds")

get_gene_scores <- function(lib){
  
  qcdf <- fread(paste0("../data/singlecell_sumstats/", exp,"_v12-mtMask_singlecell.csv.gz"), header = TRUE, sep = ",") %>% 
    data.frame() %>% filter(cell_id != "None")
  cells <- as.character(qcdf$barcode)
  fragment_file <- paste0("../../../mtscATACpaper_large_data_files/source/cellranger_output/", exp, "_v12-mtMask_fragments.tsv.gz")

  # create a gene by cell matrix
  gene.activities <- FeatureMatrix(
    fragments = fragment_file,
    features = genebodyandpromoter.coords,
    cells = cells,
    chunk = 10
  )
  
  # convert rownames from chromsomal coordinates into gene names
  gene.key <- genebodyandpromoter.coords$gene_name
  names(gene.key) <- GRangesToString(grange = genebodyandpromoter.coords)
  rownames(gene.activities) <- gene.key[rownames(gene.activities)]
  saveRDS(gene.activities, file =  paste0("../../../mtscATACpaper_large_data_files/intermediate/",lib,".gene_activities.rds"))
}
get_gene_scores("CD34_H8")
get_gene_scores("CD34_G10")
get_gene_scores("PBMC_H10")
get_gene_scores("PBMC_H9")


