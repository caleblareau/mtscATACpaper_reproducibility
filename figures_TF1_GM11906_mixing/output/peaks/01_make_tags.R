library(data.table)
library(SummarizedExperiment)
library(Matrix)
library(rtracklayer)

source("../../code/02a_assign_pull_function.R")

get_filtered_fragments <- function(frags_path, df){
  
  tf1_bcs <- df %>% filter(assign == "TF1") %>% pull(cell_id)
  merrf_bcs <- df %>% filter(assign == "GM11906") %>% pull(cell_id)
  
  frags <- fread(paste0(dir, "/outs/fragments.tsv.gz")) %>% filter(V4 %in% c(tf1_bcs, merrf_bcs))
  bc_vec <- c(rep("TF1", length(tf1_bcs)), rep("MERRF", length(merrf_bcs))); names(bc_vec) <- c(tf1_bcs, merrf_bcs)
  frags$bc <- bc_vec[as.character(frags$V4)]
  return(frags)
}

# Import fragments
fixed <- rbind(get_filtered_fragments("../../../../mtscATACpaper_large_data_files/source/cellranger_output/Mix_Fix_1hr_fragments.tsv.gz",
                                      data.frame(fread("../data1_meta.tsv"))), 
               get_filtered_fragments("../../../../mtscATACpaper_large_data_files/source/cellranger_output/Mix_Fix_6hr_fragments.tsv.gz",
                                      data.frame(fread("../data6_meta.tsv"))))
regular <- get_filtered_fragments("../../../../mtscATACpaper_large_data_files/source/cellranger_output/Mix_1_fragments.tsv.gz",
                                  data.frame(fread("../data_mix1_meta.tsv")))

export_bed <- function(what, what2, input){
  gr <- data.frame(input %>% filter(bc == what))[,c(1,2,3)]
  write.table(gr, file = paste0("", what2, "_", what, ".bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

export_bed("TF1", "fixed", fixed)
export_bed("TF1", "regular", regular)
export_bed("MERRF", "fixed", fixed)
export_bed("MERRF", "regular", regular)

# Now call peaks with MACS2





