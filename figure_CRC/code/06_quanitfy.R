library(data.table)
library(dplyr)

mc <- fread("../data/mgatk_complexity/CRC_v12-mtMask_mgatk_mitoComplexity.tsv.gz")

mc %>% dplyr::filter(barcode %in% colnames(crc)) %>% summarize( median(as.numeric(as.character(PERCENT_DUPLICATION)), NA.rm = TRUE))
