library(Matrix)
library(SummarizedExperiment)

# Function to pull the per-cell count vector of the position / letter for an mtDNA mutation
extractme <- function(position, letter){
  (assays(SE)[[paste0(letter, "_counts_fw")]][position, ] + assays(SE)[[paste0(letter, "_counts_rev")]][position, ])
}

# Function to extract the mutation and coverage for selected positions and ref / alt alleles
e2 <- function(position, letter1, letter2){
  df <- data.frame(
    m1 = extractme(position, letter1),
    m2 = extractme(position, letter2)
  )
  df$het <- round(df$m2 / (df$m2 + df$m1 + 0.0001), 3)
  df$cov <- df$m1 + df$m2
  colnames(df) <- c(paste0("m", as.character(position), "_", letter1), paste0("m", as.character(position), "_", letter2),
                    paste0("het", as.character(position), letter1, "_", letter2),
                    paste0("cov", as.character(position)))
  return(df)
}
