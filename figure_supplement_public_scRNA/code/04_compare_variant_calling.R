library(precrec)
library(Matrix)
library(SummarizedExperiment)
library(data.table)

"%ni%" <- Negate("%in%")

# Import data from 03 script
se_d1 <- readRDS(paste0("../../../mtscATACpaper_large_data_files/source/other/Donor1_smartseq2_SE_af.rds"))
se_d2 <- readRDS(paste0("../../../mtscATACpaper_large_data_files/source/other/Donor2_smartseq2_SE_af.rds"))
af_d1 <- assays(se_d1)[["allele_frequency"]]
af_d2 <- assays(se_d2)[["allele_frequency"]]

# Get homoplasmic variants to filter from other tools (freebays, bcftools)
homoplasmic1 <- data.frame(rowData(se_d1)) %>% filter(mean > 0.98) %>% pull(variant)
homoplasmic2 <- data.frame(rowData(se_d2)) %>% filter(mean > 0.98) %>% pull(variant)

# Mini function to import variants from other tools via TSV and remove any donor homoplasmic variants
process_vcf_tsv <- function(dt, homoplasmic, qual_val = 100){
  
  dt_filt <- dt %>% filter(qual >= qual_val & nchar(ref) == 1 & nchar(alt) == 1) # focus on SNVs
  # already manually verified that none are multi alt alleles
  var <- paste0(as.character(dt_filt[["start"]]), as.character(dt_filt[["ref"]]), ">",  as.character(dt_filt[["alt"]]))
  var[var %ni% homoplasmic]
}

# Plain text column names
ptcn <- c("chr", "start", "ref", "alt", "qual")
# Import processed variant calls from other tools
d1_bcf <- process_vcf_tsv(fread("../data/other_tools_variantcalls/Donor1_bcftools.tsv", skip = 1, header = FALSE, col.names = ptcn), homoplasmic1)
d1_fb <- process_vcf_tsv(fread("../data/other_tools_variantcalls/Donor1_freebayes.tsv", skip = 1, header = FALSE, col.names = ptcn), homoplasmic1)
d2_bcf <- process_vcf_tsv(fread("../data/other_tools_variantcalls/Donor2_bcftools.tsv", skip = 1, header = FALSE, col.names = ptcn), homoplasmic2)
d2_fb <- process_vcf_tsv(fread("../data/other_tools_variantcalls/Donor2_freebayes.tsv", skip = 1, header = FALSE, col.names = ptcn), homoplasmic2)

# Import variants from mgatk and Ludwig et al 2019
d1_sup <- readRDS("../output/Donor1_supervised_variants.rds") %>% pull(variant) %>% as.character()
d2_sup <- readRDS("../output/Donor2_supervised_variants.rds") %>% pull(variant) %>% as.character()
d1_mgatk <- readRDS("../output/Donor1_all_mgatk_variants.rds") %>% pull(variant) %>% as.character()
d2_mgatk <- readRDS("../output/Donor2_all_mgatk_variants.rds") %>% pull(variant) %>% as.character()


# Function to compute cell-cell mtDNA distance as well as whether or not they are in the same colony
compute_mito_dist <- function(variants, af){
  dm <- data.matrix(af[variants,])
  dm[dm < 0.01] <- 0
  dm[dm > 0.20] <- 0.20
  metric <- as.numeric(dist(t(dm)))
  colony <- str_split_fixed(colnames(af), "_", 3)[,2]
  colonymat <-   model.matrix(~0+colony)
  data.frame(same_colony = (-1*as.numeric(dist(colonymat, method = "binary"))) +1, 
             similarity = 1 + -1*metric)
}

compute_mito_dist_cor <- function(variants, af){
  
  # Subset variants
  dm <- data.matrix(af[variants,])
  
  # compute the cell-cell correlation
  cor_mat <- cor(dm)
  cor_mat[upper.tri(cor_mat, diag = TRUE)] <- NA
  cor_vec <- as.numeric(cor_mat)
  cor_vec <- cor_vec[!is.na(cor_vec)]
  
  # Append whether or not the cell was from the same colony
  colony <- str_split_fixed(colnames(af), "_", 3)[,2]
  colonymat <-   model.matrix(~0+colony)
  odf <- data.frame(same_colony = (-1*as.numeric(dist(colonymat, method = "binary"))) +1, 
                    similarity = cor_vec)
  odf
}


#-----
# Process for donor 1
#-----

d1_outcome <- compute_mito_dist_cor(d1_mgatk, af_d1)
scores_d1 <- join_scores(compute_mito_dist_cor(d1_bcf, af_d1)[,2], compute_mito_dist_cor(d1_fb, af_d1)[,2], 
                         compute_mito_dist_cor(d1_mgatk, af_d1)[,2], compute_mito_dist_cor(d1_sup, af_d1)[,2])

mmpr1 <- mmdata(scores_d1, d1_outcome$same_colony, modnames = c("bcftools", "FreeBayes", "mgatk", "Supervised"))
mscurves1 <- evalmod(mmpr1)

dfo1 <- auc(mscurves1)
dfo1$sample <- "Donor1"
dfo1

#-----
# Process for donor 2
#-----

d2_outcome <- compute_mito_dist_cor(d2_mgatk, af_d2)
scores_d2 <- join_scores(compute_mito_dist_cor(d2_bcf, af_d2)[,2], compute_mito_dist_cor(d2_fb, af_d2)[,2], 
                         compute_mito_dist_cor(d2_mgatk, af_d2)[,2], compute_mito_dist_cor(d2_sup, af_d2)[,2])

mmpr2 <- mmdata(scores_d2, d2_outcome$same_colony, modnames = c("bcftools", "FreeBayes", "mgatk", "Supervised"))
mscurves2 <- evalmod(mmpr2)

dfo2 <- auc(mscurves2)
dfo2$sample <- "Donor2"
dfo2

#-----
# Make output plots
#-----

p1 <- autoplot(mscurves1, "ROC") +
  pretty_plot(fontsize = 8) + L_border() +
  ggtitle(label = "Donor 1 colonies") +
  scale_color_manual(values = c("dodgerblue3", "purple3", "firebrick", "orange2")) +
  theme(legend.position = "none")

p2 <- autoplot(mscurves2, "ROC") +
  pretty_plot(fontsize = 8) + L_border() +
  ggtitle(label = "Donor 2 colonies") +
  scale_color_manual(values = c("dodgerblue3", "purple3", "firebrick", "orange2")) +
  theme(legend.position = "none")

p3 <- ggplot(rbind(dfo1, dfo2) %>% filter(curvetypes == "ROC"), aes(x = sample, y = aucs, fill = modnames, group = modnames)) +
  geom_bar(stat = "identity", color = "black", width=.8, position = "dodge") +
  scale_fill_manual(values = c("dodgerblue3", "purple3", "firebrick", "orange2")) +
  pretty_plot(fontsize = 8) + L_border() + scale_y_continuous(expand = c(0,0), limits = c(0,1)) +
  labs(x = "", y = "AUROC", fill = "Variants") + ggtitle(label = ".")

cowplot::ggsave(cowplot::plot_grid(p1, p2, p3, nrow = 1, rel_widths = c(1, 1, 1.5)), file = "out_plots/ROC_comparison_variant_callers.pdf", width = 7.5, height = 2)
