library(data.table)
library(dplyr)

process_mito_picard <- function(directory){
  barcodes <- gsub(".rmdups.log", "", list.files(paste0(directory, "/logs/rmdupslogs")))
  
  lapply(barcodes, function(barcode){
    dt <- fread(paste0(directory,"/logs/rmdupslogs/",barcode,".rmdups.log"))[,c("PERCENT_DUPLICATION","ESTIMATED_LIBRARY_SIZE")]
    data.frame(barcode, dt)
  }) %>% rbindlist() %>% data.frame() -> df
  write.table(df, file = paste0(directory, "_mitoComplexity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  directory
}


