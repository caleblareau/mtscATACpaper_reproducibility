library(Matrix)
library(GenomicRanges)
library(data.table)
library(dplyr)
set.seed(1)

# Perform Chi-squared test for observed and permuted chromatin accessibility values
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

process_patient_X2 <- function(patient){
  dir <- paste0("../../../mtscATACpaper_large_data_files/source/cellranger_output/CLL_",patient,"_CD19pos_peak_bc_matrix/")
  peaks <- fread(paste0(dir, "peaks.bed.gz"), header = FALSE, col.names = c("chr", "start", "end"))
  peaks_gr <- makeGRangesFromDataFrame(peaks)
  barcodes <- fread(paste0(dir, "barcodes.tsv.gz"), header = FALSE)[[1]]
  mtx <- fread(paste0(dir, "/matrix.mtx.gz")[[1]], skip = 3, header = FALSE)
  
  # Assemble binary matrix
  mat <- Matrix::sparseMatrix(i = c(mtx[[1]],length(peaks)), 
                              j = c(mtx[[2]],length(barcodes)),
                              x = c(as.numeric(mtx[[3]]) > 0, 0))
  colnames(mat) <- barcodes
  
  # Import cluster annotation from mtDNA
  cluster_df <- readRDS(paste0("../output/",patient,"_clone_definition.rds"))
  boo <- Matrix::rowSums(mat) >= 50
  mat <- mat[boo,as.character(cluster_df$cell_id)]
  clusters <- as.character(cluster_df$cluster_id)
  peaks <- peaks[boo, ]
  
  # Reformat for faster subsetting
  sm <- data.table(summary(mat))
  
  # Look at all peaks
  set.seed(1)
  x2df <- lapply(1:dim(mat)[1], function(idx){ 
    if((idx %% 1000) == 0){
      print(idx)
    }
    Y = 1:length(clusters) %in% sm[i == idx][["j"]]
    data.frame(peaks[idx,], getX2stats(Y, clusters))
    
  }) %>% rbindlist() %>% data.frame()
  
  saveRDS(x2df, file = paste0("../output/",patient,"_X2.rds"))
}

#process_patient_X2("PT1")
process_patient_X2("PT2")

