library(Matrix)
library(GenomicRanges)
library(data.table)
library(dplyr)

getX2stats <- function(y, cid){
  obs <- chisq.test(cid, y)
  perm <- chisq.test(cid, sample(y))
  data.frame(
    obs_p = obs$p.value,
    obs_x2 = obs$statistic,
    perm_p = perm$p.value,
    perm_x2 = perm$statistic
  )
}

dir <- "../CLL_CD19pos_CR-mtMask/"
peaks <- fread(paste0(dir, "/outs/filtered_peak_bc_matrix/peaks.bed"), header = FALSE, col.names = c("chr", "start", "end"))
peaks_gr <- makeGRangesFromDataFrame(peaks)
barcodes <- fread(paste0(dir, "/outs/filtered_peak_bc_matrix/barcodes.tsv"), header = FALSE)[[1]]
mtx <- fread(paste0(dir, "/outs/filtered_peak_bc_matrix/matrix.mtx")[[1]], skip = 3, header = FALSE)

# Assemble binary matrix
mat <- Matrix::sparseMatrix(i = c(mtx[[1]],length(peaks)), 
                            j = c(mtx[[2]],length(barcodes)),
                            x = c(as.numeric(mtx[[3]]) > 0, 0))
colnames(mat) <- barcodes

# Import cluster annotation from mtDNA
cluster_df <- readRDS("output/PT1_clone_definition.rds")
boo <- Matrix::rowSums(mat) >= 50
mat <- mat[boo,as.character(cluster_df$cell_id)]
clusters <- as.character(cluster_df$cluster_id)
peaks <- peaks[boo, ]

# Reformat for faster subsetting
sm <- data.table(summary(mat))


set.seed(1)
x2df <- lapply(1:dim(mat)[1], function(idx){ 
  if((idx %% 1000) == 0){
    print(idx)
  }
  Y = 1:length(clusters) %in% sm[i == idx][["j"]]
  data.frame(peaks[idx,], getX2stats(Y, clusters))
  
}) %>% rbindlist() %>% data.frame()

saveRDS(x2df, file = "output/PT1_X2.rds")
