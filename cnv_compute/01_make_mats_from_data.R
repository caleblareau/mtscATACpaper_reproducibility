library(data.table)
library(dplyr)
source("00_cnv_matrix_make.R")

# Create baseline CNV matrix from public, healthy data
mat <- make_basic_mat("../../mtscATACpaper_large_data_files/source/other/public_atac_data_10x/atac_v1_pbmc_10k_fragments.tsv.gz",
                      fread("../../mtscATACpaper_large_data_files/source/other/public_atac_data_10x/atac_v1_pbmc_10k_barcodes.tsv", header = FALSE)[[1]])

saveRDS(mat, file = "output/atac_public_pbmcs_cnv.rds")

# Quick helper function to grab passed barcodes from cellranger
pull_barcodes_sc_file <- function(sc_file){
  data.frame(fread(sc_file)) %>% filter(cell_id != "None") %>% pull(barcode)
}

# Process CRC data
CRC_barcodes <- pull_barcodes_sc_file("../figure_CRC/data/singlecell_sumstats/CRC_singlecell.csv.gz")
mat <- make_basic_mat("../../mtscATACpaper_large_data_files/source/cellranger_output/CRC_v12-mtMask_fragments.tsv.gz",
                      CRC_barcodes)
saveRDS(mat, file = "output/CRC_tumor_cnv.rds")

CLLPT1_CD19pos_barcodes <- pull_barcodes_sc_file("../figure_CLL/data/singlecell_sumstats/CLL_PT1_CD19pos_v12-mtMask_singlecell.csv.gz")
mat <- make_basic_mat("../../mtscATACpaper_large_data_files/source/cellranger_output/CLL_PT1_CD19pos_v12-mtMask_fragments.tsv.gz",
                      CLLPT1_CD19pos_barcodes)
saveRDS(mat, file = "output/CLL_PT1_CD19pos_cnv.rds")

