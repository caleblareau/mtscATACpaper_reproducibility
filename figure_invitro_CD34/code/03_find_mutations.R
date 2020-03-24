library(SummarizedExperiment)
library(Matrix)
library(dplyr)
library(BuenColors)
library(data.table)
source("../../global_functions/variant_calling.R")
load("../output/trajectory_inferences.18march2020.rda")

"%ni%" <- Negate("%in%")
perc.rank <- function(x) trunc(rank(x))/length(x)

# Import SEs
import_SE <- function(str){
  SE <- readRDS(paste0("../../../mtscATACpaper_large_data_files/source/mgatk_output/",str,"_v12-mtMask_mgatk.rds"))
  colnames(SE) <- paste0(str, "-", as.character(colnames(SE)))
  SE <- SE[,rownames(plot_df)[which(plot_df$library == str)]]
  return(SE)
}
SE500 <- SummarizedExperiment::cbind(import_SE("CD34_500_Day08"), import_SE("CD34_500_Day14"))
SE800 <- SummarizedExperiment::cbind(import_SE("CD34_800_Day08"), import_SE("CD34_800_Day14"), import_SE("CD34_800_Day20"))

process_clonal_experiment<- function(SE, countc){
  
  # Call mutations
  mut_se <- call_mutations_mgatk(SE, stabilize_variance = FALSE)
  pp_all <- data.frame(rowData(mut_se))
  
  transition <- c("C>T", "G>A", "A>G", "T>C")
  pp_all$transition <- pp_all$nucleotide %in% transition
  
  p1 <- ggplot(shuf(pp_all) %>% filter(n_cells_conf_detected >= 5), aes(x = strand_correlation, y = log10(vmr), color = transition)) +
    geom_point(size = 0.02) + scale_color_manual(values = c("black", "firebrick")) +
    labs(color = "HQ", x = "Strand concordance", y = "log VMR") +
    pretty_plot(fontsize = 7) + L_border() +
    theme(legend.position = "bottom") + geom_vline(xintercept = 0.65, linetype = 2) + 
    geom_hline(yintercept = -2, linetype = 2) + theme(legend.position = "none")
  cowplot::ggsave2(p1, file = paste0("../plots/call_vars_",countc,"in.pdf"), width = 1.2, height = 1.2)
  
  pp <- pp_all %>% filter(n_cells_conf_detected >= 5 & log10(vmr) > -2 & strand_correlation >= 0.65)
  sum(pp$transition)/length(pp$nucleotide)

  mut_se2 <- mut_se[pp$variant,]
  saveRDS(mut_se2, file = paste0("../output/filtered_mitoSE_CD34-",countc,".rds"))
  dim(pp)
}
process_clonal_experiment(SE500, "500")
process_clonal_experiment(SE800, "800")
