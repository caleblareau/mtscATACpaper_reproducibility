library(data.table)
library(BuenColors)
library(dplyr)
library(ggbeeswarm)
library(SummarizedExperiment)

# count fraction of cells retained by mito filter
rbind(fread("../output/barcode_qc/CD34_G10.barcode_qc.tsv"),
      fread("../output/barcode_qc/CD34_H8.barcode_qc.tsv")) %>% filter(FRIP > 0.25 & depth > 1000) %>% pull(keep) %>% mean

rbind(fread("../output/barcode_qc/PBMC_H10.barcode_qc.tsv"),
      fread("../output/barcode_qc/PBMC_H9.barcode_qc.tsv")) %>% filter(FRIP > 0.60 & depth > 1000) %>% pull(keep) %>% mean


# Pull meta data per cell  from cellranger output
import_sc <- function(lib){
  dt <- fread(paste0("../data/singlecell_sumstats/",lib,"_v12-mtMask_singlecell.csv.gz"))
  dt %>% filter(cell_id != "None") %>% mutate(pct_mito = mitochondrial/total*100) %>%
    mutate(barcode = paste0(lib, "_", barcode)) %>% data.frame()
}

# Pull meta data per cell from mgatk output via Picard
import_mk <- function(lib){
  dt <- fread(paste0("../data/mgatk_complexity/",lib,"_v12-mtMask_mgatk_mitoComplexity.tsv.gz"))
  dt %>%  mutate(barcode = paste0(lib, "_", barcode)) %>% data.frame()
}

# Import for the PBMC libraries
single_cell_metrics <- rbind(import_sc("PBMC_H10"),import_sc("PBMC_H9"))
mgatk_complexity_metrics <- rbind(import_mk("PBMC_H10"),import_mk("PBMC_H9"))
sdf <- data.frame(colData(readRDS("../output/filteredpbmcs_mgatk_calls.rds"))); sdf$barcode <- rownames(sdf)
mdf_pbmc <- merge(merge(sdf, single_cell_metrics, by = "barcode"), mgatk_complexity_metrics, by = "barcode")

# import for the CD34 libraries
single_cell_metrics <- rbind(import_sc("CD34_H8"),import_sc("CD34_G10"))
mgatk_complexity_metrics <- rbind(import_mk("CD34_H8"),import_mk("CD34_G10"))
sdf <- data.frame(colData(readRDS("../output/filteredCD34_mgatk_calls.rds"))); sdf$barcode <- rownames(sdf)
mdf_cd34 <- merge(merge(sdf, single_cell_metrics, by = "barcode"), mgatk_complexity_metrics, by = "barcode")

# Combine for plotting
mdf_pbmc$what <- "PBMC"
mdf_cd34$what <- "CD34"
mdf_plot <- rbind(mdf_cd34, mdf_pbmc)
mdf_plot %>% group_by(what) %>% summarize(mean(passed_filters), median(as.numeric(as.character(PERCENT_DUPLICATION)), NA.rm = TRUE))


# Visualize the percent duplication of reads
p1 <- ggplot(mdf_plot, aes(x = what, y = as.numeric(as.character(PERCENT_DUPLICATION))*100)) +
  geom_boxplot(outlier.shape = NA) + scale_y_continuous(limits = c(10, 60)) +
  pretty_plot(fontsize = 7) + L_border() +
  labs(x = "", y = "% duplicates in mtDNA")

# Visualize the fraction of mtDNA reads per celltype
p2 <- ggplot(mdf_plot, aes(x = what, y = depth)) +
  geom_boxplot(outlier.shape = NA) + scale_y_continuous(limits = c(20, 150)) +
  pretty_plot(fontsize = 7) + L_border() +
  labs(x = "", y = "Mean mtDNA coverage")

p3 <- ggplot(mdf_plot, aes(x = what, y = pct_mito)) +
  geom_boxplot(outlier.shape = NA) + scale_y_continuous(limits = c(10, 60)) +
  pretty_plot(fontsize = 7) + L_border() +
  labs(x = "", y = "% mitochondrial reads")

cowplot::ggsave2(cowplot::plot_grid(p1, p2, p3, nrow = 1), width = 2.4, height = 1.7, filename = "../plots/CD34_PBMC_compare_QC.pdf")
