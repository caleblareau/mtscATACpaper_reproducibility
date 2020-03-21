library(data.table)
library(Matrix)

# From manual curation of previous heatmap
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
    tf1_mean = round(colMeans(mmat[tf1_vars,]), 3),
    gm11906_mean = round(colMeans(mmat[gm11906_vars,]), 3),
    merrf_af  = round(mmat[merrf,], 3),
    merrf_cov = round(cov[8344,], 0),
    mean_cov = round(colMeans(cov), 1)
  )
  df$assign <- ifelse(
    df$gm11906_mean > 0.95, "GM11906", 
    ifelse(df$tf1_mean > 0.95, "TF1", 
           ifelse(df$mean_cov < 5, "Low_coverage", "Collision")))
  df
  
}