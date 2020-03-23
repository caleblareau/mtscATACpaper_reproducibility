library(SummarizedExperiment)
library(Matrix)
library(dplyr)
library(BuenColors)
library(ComplexHeatmap)
library(circlize)
library(data.table)
library(ggrepel)
source("../../global_functions/variant_calling.R")

"%ni%" <- Negate("%in%")
perc.rank <- function(x) trunc(rank(x))/length(x)

# Import mgatk object
import_mgatk <- function(library){
  df <- fread(paste0("../output/barcode_qc/",library,".barcode_qc.tsv")) %>%
    filter(keep)
  se <- readRDS(paste0("../../../mtscATACpaper_large_data_files/source/mgatk_output/",library,"_v12-mtMask_mgatk.rds"))
  se_filt <- se[,as.character(df$sample)]
  colnames(se_filt) <- paste0(library, "_", colnames(se_filt))
  se_filt
}

#--------------------
# Call CD34 mutations
#--------------------
SE_CD34 <- cbind(import_mgatk("CD34_G10"), import_mgatk("CD34_H8"))

# Call variants
mut_se_CD34 <- call_mutations_mgatk(SE_CD34)
misc_df_CD34 <- data.frame(rowData(mut_se_CD34))
vars_cd34 <- misc_df_CD34 %>%  filter(n_cells_detected >= 5 & strand_correlation > 0.65 & log10(vmr) > -2 & mean_coverage >= 20) %>% pull(variant)

# Mark transitions
transition <- c("C>T", "G>A", "A>G", "T>C")
misc_df_CD34$transition <- misc_df_CD34$nucleotide %in% transition

pC <- ggplot(misc_df_CD34 %>%  filter(n_cells_detected >= 5 & mean_coverage >= 20), aes(x = strand_correlation, y = log10(vmr), color = transition)) +
  geom_point(size = 0.4) + scale_color_manual(values = c("black", "dodgerblue")) +
  labs(color = "HQ", x = "Strand concordance", y = "log VMR") +
  pretty_plot(fontsize = 8) + L_border() +
  theme(legend.position = "bottom") + 
  geom_vline(xintercept = 0.65, linetype = 2) +
  geom_hline(yintercept = -2, linetype = 2) + theme(legend.position = "none")
cowplot::ggsave2(pC, file = "../plots/invivoCD34_called_variants.pdf", width = 1.7, height = 1.7)


#--------------------
# Call PBMC mutations
#--------------------
SE_PBMC <- cbind(import_mgatk("PBMC_H9"), import_mgatk("PBMC_H10"))

# Call variants
mut_se_pbmc <- call_mutations_mgatk(SE_PBMC)
misc_df_pbmc <- data.frame(rowData(mut_se_pbmc))
vars_pbmc <- misc_df_pbmc %>%  filter(n_cells_detected >= 5 & strand_correlation > 0.65 & log10(vmr) > -2 & mean_coverage >= 20) %>% pull(variant)

# Mark transitions
transition <- c("C>T", "G>A", "A>G", "T>C")
misc_df_pbmc$transition <- misc_df_pbmc$nucleotide %in% transition

pP <- ggplot(misc_df_pbmc %>%  filter(n_cells_detected >= 5 & mean_coverage >= 20), aes(x = strand_correlation, y = log10(vmr), color = transition)) +
  geom_point(size = 0.4) + scale_color_manual(values = c("black", "dodgerblue")) +
  labs(color = "HQ", x = "Strand concordance", y = "log VMR") +
  pretty_plot(fontsize = 8) + L_border() +
  theme(legend.position = "bottom") + 
  geom_vline(xintercept = 0.65, linetype = 2) +
  geom_hline(yintercept = -2, linetype = 2) + theme(legend.position = "none")
cowplot::ggsave2(pP, file = "../plots/invivoPBMCs_called_variants.pdf", width = 1.7, height = 1.7)

# Export the variants called in either culture
length(vars_pbmc)
length(vars_cd34)
vars_both <- unique(c(vars_pbmc, vars_cd34))

mscf <- mut_se_CD34[vars_both,]
mspf <- mut_se_pbmc[vars_both,]

# Annotate where the variants came from 
rowData(mscf)$called_in_CD34 <- rownames(mscf) %in% vars_cd34
rowData(mscf)$called_in_PBMC <- rownames(mscf) %in% vars_pbmc
rowData(mspf)$called_in_CD34 <- rownames(mspf) %in% vars_cd34
rowData(mspf)$called_in_PBMC <- rownames(mspf) %in% vars_pbmc

saveRDS(mscf, file = "../output/filteredCD34_mgatk_calls.rds")
saveRDS(mspf, file = "../output/filteredpbmcs_mgatk_calls.rds")

mean_df <- data.frame(
  CD34 = rowData(mscf)$mean,
  PBMC = rowData(mspf)$mean, 
  variant = rownames(mscf),
  v2 = rownames(mspf)
) %>% mutate(log10_FC = log10(CD34/PBMC))  %>% arrange(desc(log10_FC))

ggplot(mean_df, aes(x = CD34, y = PBMC, color = variant %in% c("3244G>A", "11318T>G", "174C>T"))) +
  geom_point() + scale_x_log10(limits = c(1e-04, 0.02)) + scale_y_log10(limits = c(1e-04, 0.02)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) + 
  pretty_plot() + L_border()
