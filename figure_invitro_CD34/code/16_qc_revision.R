library(data.table)
library(BuenColors)
library(dplyr)
library(ggbeeswarm)
library(SummarizedExperiment)

# count fraction of cells retained by mito filter
rbind(fread("../output/barcode_qc/CLL_PT1_CD19pos.barcode_qc.tsv"),
      fread("../output/barcode_qc/CLL_PT2_CD19pos.barcode_qc.tsv")) %>% filter(FRIP > 0.4 & depth > 1000) %>% mutate(keep = mtDNAcoverage > 20) %>% pull(keep) %>% mean


#  cell from mgatk output via Picard
import_mk <- function(lib){
  dt <- fread(paste0("../data/mgatk_complexity/",lib,"_v12-mtMask_mgatk_mitoComplexity.tsv.gz"))
  dt %>%  mutate(barcode = paste0(lib, "-", barcode)) %>% data.frame()
}

# Import for the PBMC libraries
mgatk_complexity_metrics <- rbind(import_mk("CLL_PT1_CD19pos"),import_mk("CLL_PT2_CD19pos"))

merge <- merge(sdf, mgatk_complexity_metrics, by = "barcode")

