library(data.table)
library(BuenColors)
library(dplyr)
library(ggbeeswarm)

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

# Import the HD cell type labels
sdf <- readRDS("../output/PBMCatac_SignacSeurat_labelTransfer.rds")
sdf$barcode <- rownames(sdf)
mdf <- merge(merge(sdf, single_cell_metrics, by = "barcode"), mgatk_complexity_metrics, by = "barcode")

# Filter out annotations that have less ten cells
mdf_plot  <- mdf %>% group_by(predicted.id) %>% mutate(count = n()) %>% dplyr::filter(count > 10)
mdf_plot$predicted.id <- gsub("cyte", "", gsub("cell", "", as.character(mdf_plot$predicted.id)), "cyte")
mdf_plot$predicted.id <- factor(mdf_plot$predicted.id,
                                levels = c("Memory_CD4_T","Naive_CD4_T","Treg", "GammaDelta_T","Cytotoxic_CD8_T","NK_",
                                           "activated_B", "Memory_B","CD14_mono", "CD16_mono","Dendritic_", "pDC"))

# Visualize the percent duplication of reads
p1 <- ggplot(mdf_plot, aes(x = predicted.id, y = PERCENT_DUPLICATION*100)) +
  geom_boxplot(outlier.shape = NA) + scale_y_continuous(limits = c(10, 60)) +
  pretty_plot(fontsize = 7) + L_border() +
  labs(x = "", y = "% duplicates in mtDNA")

# Visualize the fraction of mtDNA reads per celltype
p2 <- ggplot(mdf_plot, aes(x = predicted.id, y = pct_mito)) +
  geom_boxplot(outlier.shape = NA) + scale_y_continuous(limits = c(10, 60)) +
  pretty_plot(fontsize = 7) + L_border() +
  labs(x = "", y = "% mitochondrial reads")

#-----
# Process the RNA-seq data from external resource as replication
#-----
mito_genes <- fread("../data/haemosphere_rnaseq/Homo_sapiens.gene_info.gz") %>% 
  filter(chromosome == "MT") %>% pull(GeneID)
mito_genes_ens <- fread("../data/haemosphere_rnaseq/AllGenes.txt.gz") %>% filter(Species == "HomoSapiens") %>%
  filter(EntrezId %in% mito_genes) %>% pull(EnsemblId)
tpm <- fread("../data/haemosphere_rnaseq/Haemopedia-Human-RNASeq_tpm.txt.gz") %>%
  filter(genenames %in% mito_genes_ens)

# Add mtDNA abundances to 
total_mito_gene_tpm <- colSums(tpm[,-1])
d <- total_mito_gene_tpm[order(names(total_mito_gene_tpm))]
dd <-fread("../data/haemosphere_rnaseq/Haemopedia-Human-RNASeq_samples.txt.gz") %>% arrange(sampleId) %>%
  mutate(tpm_pct = d/1000000*100)

# Collapse some meta data
dd$cta <- case_when(
  dd$cellTypeAbbreviation %in% c("CD4T", "CD8T") ~ "Tcell",
  dd$cellTypeAbbreviation %in% c("MemB", "NveB") ~ "Bcell",
  dd$cellTypeAbbreviation %in% c("MonoNonClassic ", "Mono") ~ "Mono",
  dd$cellTypeAbbreviation %in% c("myDC", "myDC123") ~ "Dendritic",
  dd$cellTypeAbbreviation %in% c("NK") ~ "NK",
  dd$cellTypeAbbreviation %in% c("pDC") ~ "pDC",
  TRUE ~ "other"
)

dd$cta <- factor(dd$cta,
                 levels = c("Tcell", "NK", "Bcell", "Mono", "Dendritic", "pDC", "other"))


p3 <- ggplot(dd %>% filter(!(cta %in% c("other"))), aes(x = cta, y = tpm_pct)) +
  geom_boxplot(outlier.shape =  NA)+ pretty_plot(fontsize = 7) + L_border() +
  labs(x = "", y = "% mitochondrial counts")


cowplot::ggsave2(cowplot::plot_grid(p1, p2, p3, nrow = 1, rel_widths =c(1,1,0.5)),
                 width = 7.5, height = 1.5, file = "../plots/populations_duplicate_rate.pdf")

