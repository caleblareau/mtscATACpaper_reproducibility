library(data.table)
library(Matrix)
library(BuenColors)
library(irlba)
library(umap)
library(dplyr)
set.seed(1)

# Quick script to subset genes into 24 bins for module scoring
# The idea here is to use an external healthy control
existing_genes <- make.unique(fread("../../scRNAseq_data/gene_expression/CLL1_CD19neg_ge/features.tsv.gz", header = FALSE)[[2]])

# Imoprt data from the 10x v3 5k platform
dt <- fread("zcat < data/5k_v3/matrix.mtx.gz", skip = 3); barcodes <- fread("zcat < data/5k_v3/barcodes.tsv.gz", header = FALSE)[[1]]
genes <- make.unique(fread("data/5k_v3/features.tsv.gz", header = FALSE)[[2]])
mat <- sparseMatrix(
  i = c(dt[["V1"]],length(genes)),
  j = c(dt[["V2"]],1),
  x = c(dt[["V3"]], 0)
)
rownames(mat) <- genes; colnames(mat) <- barcodes
mat <- mat[genes %in% existing_genes,]
count_norm <- log1p(t(t(mat)/Matrix::colSums(mat)) * 10000)

# Find expressed genes... 12,713
cpm_gene <- rowSums(mat)/sum(rowSums(mat)) *1000000
sorted_gene_vec <- sort(cpm_gene[cpm_gene >= 1]) 
nbin <- 24
bin_id <- reshape2::melt(split(sorted_gene_vec, ceiling(seq_along(sorted_gene_vec)/ceiling(length(sorted_gene_vec)/nbin))))
bin_id$gene <- names(sorted_gene_vec)

# Export the bin data frame in order to score the CLL samples
saveRDS(bin_id, file = "bin_df_PBMCs_healthy.rds")


#---
# The rest of the code is not used for the paper but serves as a sanity check that the module scoring worked
#--- 

# score gene signatures
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


# Import gene sets
bcell_genes <- fread("../ICA_marker_genes.tsv") %>% filter(celltype == "Follicular_B_cell") %>% pull(marker_gene) %>% unique()
myeloid_genes <- fread("../ICA_marker_genes.tsv") %>% filter(celltype %in% c("Monocyte", "Dendritic_Cell")) %>% pull(marker_gene) %>% unique()
tcell_genes <- fread("../ICA_marker_genes.tsv") %>% filter(celltype %in% c("CD8_T_cell", "Naive_T_cell")) %>% pull(marker_gene) %>% unique()
nk_genes <- fread("../ICA_marker_genes.tsv") %>% filter(celltype == "NK_cells") %>% pull(marker_gene) %>% unique()

# Test gene scores as sanity check
bscore <- scoreCellsSignature(count_norm, bin_id, bcell_genes)
myeloid_score <- scoreCellsSignature(count_norm, bin_id, myeloid_genes)
tcell_score <- scoreCellsSignature(count_norm, bin_id, tcell_genes)
nk_score <- scoreCellsSignature(count_norm, bin_id, nk_genes)

qplot(bscore, nk_score)


