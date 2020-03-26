library(data.table)
library(SummarizedExperiment)
library(dplyr)
library(Matrix)

samples1 <- gsub("_QC.rds", "", list.files("rds_SE/")[1:6])

import_meta <- function(what, what2){
  SE <- readRDS(paste0("rds_SE/", what, "_QC.rds"))
  SE$mtDNACoverage <- fread(paste0(what2, "_CR-mtMask_mgatk/final/", what2, "_CR-mtMask_mgatk.depthTable.txt")) %>%
    dplyr::filter(V1 %in% colnames(SE)) %>% pull(V2)
  colData(SE)
}

df <- rbind(
  import_meta("Patient1_CLL_15min", "CLL_15min"), import_meta("Patient1_CLL_5min", "CLL_5min"),
  import_meta("Patient1_CLL_neg", "CLL_CD19neg"), import_meta("Patient1_CLL_pos", "CLL_CD19pos"),
  import_meta("Patient2_CLL_neg", "Patient2_CLL_CD19_neg"), import_meta("Patient2_CLL_pos", "Patient2_CLL_CD19_pos")
)

fdf <- data.frame(df) %>% dplyr::filter(mtDNACoverage >= 20)
