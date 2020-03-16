library(data.table)
library(Matrix)
library(BuenColors)
library(irlba)
library(umap)
library(dplyr)
set.seed(1)

# Pull immune cell marker genes from https://github.com/caleblareau/immune_cell_signature_genes
bcell_genes <- fread("../scRNAseq_data/ICA_marker_genes.tsv") %>% filter(celltype == "Follicular_B_cell") %>% pull(marker_gene) %>% unique()
myeloid_genes <- fread("../scRNAseq_data/ICA_marker_genes.tsv") %>% filter(celltype %in% c("Monocyte", "Dendritic_Cell")) %>% pull(marker_gene) %>% unique()
tcell_genes <- fread("../scRNAseq_data/ICA_marker_genes.tsv") %>% filter(celltype %in% c("CD8_T_cell", "Naive_T_cell")) %>% pull(marker_gene) %>% unique()
nk_genes <- fread("../scRNAseq_data/ICA_marker_genes.tsv") %>% filter(celltype == "NK_cells") %>% pull(marker_gene) %>% unique()

# Import gene bin assignments based on raw expression
bindf <- readRDS("../scRNAseq_data/healthy_control_module_scoring/bin_df_PBMCs_healthy.rds")

# score gene signatures given a target gene dataset
# This function is essentially a stripped down version of Seurat's AddModuleScore (https://rdrr.io/github/satijalab/seurat/man/AddModuleScore.html)
scoreCellsSignature <- function(count_norm, bin_df, target_genes){
  count_norm_ss <- data.matrix(count_norm)[bin_df$gene, ]
  
  # Subset target genes or not
  target_boo <- rownames(count_norm_ss) %in% target_genes
  target_mat <- count_norm_ss[target_boo,]
  offtarget_mat <- count_norm_ss[!target_boo,]
  bin_vec <- bin_df$L1; names(bin_vec) <- as.character(bin_df$gene)
  
  # Create means of bin
  bins <- unique(bin_df$L1)
  bin_mean_mat <- sapply(bins, function(bin){
    bin_genes <- bin_df[bin_df$L1 == bin, "gene"]
    colMeans(offtarget_mat[rownames(offtarget_mat) %in% bin_genes,])
  }) %>% t()
  
  bin_number <- bin_vec[rownames(target_mat)]
  signature_score <- colMeans(target_mat - bin_mean_mat[as.numeric(unname(bin_number)),])
  return(signature_score)
}

# Function to score signatures based on an RNA-seq counts matrix
project_mats_do_all <- function(mat){
  
  cpm_new <- log1p(t(t(mat)/Matrix::colSums(mat)) * 10000)
  rownames(cpm_new) <- rownames(mat); colnames(cpm_new) <- colnames(mat)
  
  pldf2 <- data.frame(
    Bcell_signature = scoreCellsSignature(cpm_new, bindf, bcell_genes),
    Tcell_signature = scoreCellsSignature(cpm_new, bindf, tcell_genes),
    NK_signature = scoreCellsSignature(cpm_new, bindf, nk_genes),
    Myeloid_signature = scoreCellsSignature(cpm_new, bindf, myeloid_genes),
    barcode = colnames(cpm_new),
    n_umis = colSums(mat), n_genes = colSums(mat > 0)
  )
  return(pldf2)
}
