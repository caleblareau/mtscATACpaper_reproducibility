library(data.table)
library(SummarizedExperiment)
library(Matrix)
library(rtracklayer)
library(dplyr)
source("02a_assign_pull_function.R")

# Function to filter fragments based on QC cells and also assign cell type for TF1/GM11906 celltypes
get_filtered_fragments <- function(mgatk_file, fragments_file){
  
  # Assign cell types
  assign_df <- assign_pull( readRDS(mgatk_file))
  tf1_bcs <- assign_df %>% filter(assign == "TF1") %>% pull(cell_id) %>% as.character()
  merrf_bcs <- assign_df %>% filter(assign == "GM11906") %>% pull(cell_id) %>% as.character()
  
  # Import fragments
  frags <- fread(fragments_file) %>% filter(V4 %in% c(tf1_bcs, merrf_bcs))
  bc_vec <- c(rep("TF1", length(tf1_bcs)), rep("MERRF", length(merrf_bcs))); names(bc_vec) <- c(tf1_bcs, merrf_bcs)
  frags$bc <- bc_vec[as.character(frags$V4)]
  return(frags)
}

# Import fragments based on mitochondria cell type assignments
fixed <- rbind(get_filtered_fragments("../../../mtscATACpaper_large_data_files/source/mgatk_output/Mix_Fix_1h_CR-mtMask_mgatk.rds",
                                      "../../../mtscATACpaper_large_data_files/source/cellranger_output/Mix_Fix_1hr_fragments.tsv.gz"), 
               get_filtered_fragments("../../../mtscATACpaper_large_data_files/source/mgatk_output/Mix_Fix_6h_CR-mtMask_mgatk.rds",
                                      "../../../mtscATACpaper_large_data_files/source/cellranger_output/Mix_Fix_6hr_fragments.tsv.gz"))
regular <- get_filtered_fragments("../../../mtscATACpaper_large_data_files/source/mgatk_output/Mix_1_CR-mtMask_mgatk.rds",
                                  "../../../mtscATACpaper_large_data_files/source/cellranger_output/Mix_1_fragments.tsv.gz")

# Create and export bigwigs for later viewing
export_bw <- function(what, what2, input){
  gr <- input %>% filter(bc == what) %>% setnames(c("chr", "start", "end", "barcode", "p", "bc")) %>%
    makeGRangesFromDataFrame()
  reads_coverage <- coverage(gr)/length(gr)*1000000
  export.bw(reads_coverage, con = paste0("../../../mtscATACpaper_large_data_files/intermediate/celltype_bigwigs/", what2, "_", what, ".bw"))
}

# Process samples
export_bw("TF1", "fixed", fixed)
export_bw("TF1", "regular", regular)
export_bw("MERRF", "fixed", fixed)
export_bw("MERRF", "regular", regular)



