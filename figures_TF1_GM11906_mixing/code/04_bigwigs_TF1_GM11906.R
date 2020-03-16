library(data.table)
library(SummarizedExperiment)
library(Matrix)
library(rtracklayer)

merrf <- "8344A>G"
tf1_vars <- c("9824T>C", "6680T>C", "8701A>G", "15301G>A", "12811T>C", "14783T>C", "4164A>G", "7853G>A", 
              "7684T>C", "10345T>C", "12705C>T", "15043G>A", "5351A>G", "10873T>C", "150C>T", "6455C>T",
              "9540T>C", "12405C>T", "199T>C", "16129G>A", "10400C>T", "8552T>C", "4071C>T", "4048G>A",
              "16298T>C", "16297T>C", "16223C>T", "5460G>A", "16189T>C")
gm11906_vars <- c("6734G>A", "3010G>A", "6554C>T", "14798T>C", "11251A>G", "15452C>A", "12612A>G", "16069C>T",
                  "16092T>C", "16126T>C", "12127G>A", "462C>T", "4216T>C", "16261C>T", "228G>A", "295C>T")


assign_pull <- function(SE){
  cov <- assays(SE)[["coverage"]]+ 0.001
  ref_allele <- as.character(rowRanges(SE)$refAllele)
  
  getMutMatrix <- function(letter){
    mat <- (assays(SE)[[paste0(letter, "_counts_fw")]] + assays(SE)[[paste0(letter, "_counts_rev")]]) / cov
    rownames(mat) <- paste0(as.character(1:dim(mat)[1]), toupper(ref_allele), ">", letter)
    return(mat[toupper(ref_allele) != letter,])
  }
  
  mmat <- rbind(getMutMatrix("A"), getMutMatrix("C"), getMutMatrix("G"), getMutMatrix("T"))
  df <- data.frame(
    cell_id = colnames(mmat),
    tf1_mean = round(Matrix::colMeans(mmat[tf1_vars,]), 3),
    gm11906_mean = round(Matrix::colMeans(mmat[gm11906_vars,]), 3)
  )
  df$assign <- ifelse(
    df$gm11906_mean > 0.95, "GM11906", 
    ifelse(df$tf1_mean > 0.95, "TF1", 
           ifelse(df$mean_cov < 5, "Low_coverage", "Collision")))
  df
  
}


get_filtered_fragments <- function(dir){
  # Assign cell types
  base_dir <- basename(dir)
  assign_df <- assign_pull( readRDS(paste0(dir, "_mgatk/final/", base_dir, "_mgatk.rds")))
  tf1_bcs <- assign_df %>% filter(assign == "TF1") %>% pull(cell_id) %>% as.character()
  merrf_bcs <- assign_df %>% filter(assign == "GM11906") %>% pull(cell_id) %>% as.character()
  
  frags <- fread(paste0(dir, "/outs/fragments.tsv.gz")) %>% filter(V4 %in% c(tf1_bcs, merrf_bcs))
  bc_vec <- c(rep("TF1", length(tf1_bcs)), rep("MERRF", length(merrf_bcs))); names(bc_vec) <- c(tf1_bcs, merrf_bcs)
  frags$bc <- bc_vec[as.character(frags$V4)]
  return(frags)
}

fixed <- rbind(get_filtered_fragments("Mix_Fix_1h_CR-mtMask"), 
               get_filtered_fragments("Mix_Fix_6h_CR-mtMask"))
regular <- get_filtered_fragments("../30October2018/Mix_1_CR-mtMask")

export_bw <- function(what, what2, input){
  gr <- input %>% filter(bc == what) %>% setnames(c("chr", "start", "end", "barcode", "p", "bc")) %>%
    makeGRangesFromDataFrame()
  reads_coverage <- coverage(gr)/length(gr)*1000000
  export.bw(reads_coverage, con = paste0("bigwigs/", what2, "_", what, ".bw"))
}

export_bw("TF1", "fixed", fixed)
export_bw("TF1", "regular", regular)
export_bw("MERRF", "fixed", fixed)
export_bw("MERRF", "regular", regular)



