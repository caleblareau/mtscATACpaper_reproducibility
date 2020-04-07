library(data.table)
library(BuenColors)
library(dplyr)
library(ggbeeswarm)
library(SummarizedExperiment)

# count fraction of cells retained by mito filter
rbind(fread("../output/barcode_qc/CD34_500_Day08.barcode_qc.tsv"),
      fread("../output/barcode_qc/CD34_500_Day14.barcode_qc.tsv"),
      fread("../output/barcode_qc/CD34_800_Day08.barcode_qc.tsv"),
      fread("../output/barcode_qc/CD34_800_Day14.barcode_qc.tsv"),
      fread("../output/barcode_qc/CD34_800_Day20.barcode_qc.tsv")) %>% filter(FRIP > 0.25 & depth > 1000) %>% mutate(keep = mtDNAcoverage > 20) %>% pull(keep) %>% mean


#  cell from mgatk output via Picard
import_mk <- function(lib){
  dt <- fread(paste0("../data/mgatk_complexity/",lib,"_v12-mtMask_mgatk_mitoComplexity.tsv.gz"))
  dt %>%  mutate(barcode = paste0(lib, "-", barcode)) %>% data.frame()
}

# Import for the PBMC libraries
mgatk_complexity_metrics <- rbind(import_mk("CD34_500_Day08"),import_mk("CD34_500_Day14"), import_mk("CD34_800_Day08"), import_mk("CD34_800_Day14"), import_mk("CD34_800_Day20"))
sdf <- rbind(data.frame((readRDS("../output/cluster-id-500.rds"))),
             data.frame((readRDS("../output/cluster-id-800.rds"))))
merge <- merge(sdf, mgatk_complexity_metrics, by = "barcode")

