library(SummarizedExperiment)
library(Matrix)
library(dplyr)
library(BuenColors)
library(ComplexHeatmap)
library(circlize)
library(data.table)
library(ggrepel)
library(ggbeeswarm)
library(stringr)
source("../../global_functions/variant_calling.R")
"%ni%" <- Negate("%in%")

# General fuction
substrRight1 <- function(x, n=1){
  substr(x, nchar(x)-n+1, nchar(x))
}

# import paired ATAC/RNA data from TF1 to identify RNA-biased variants
pseudo <- 0.0000000001
tf1_data <- readRDS("../data/meta_data/TF1_longDataFrame_ATACRNA.rds")
#rna_biased_vars <- tf1_data %>% filter(round(atac_af,2) < 0.0001 & round(rna_af, 2) > 0) %>% pull(pos) %>% as.character() %>% unique()
rna_biased_vars_10 <- tf1_data %>% mutate(ratio = (rna_af + pseudo)/(atac_af + pseudo)) %>% 
  filter(ratio > 10 )%>% # 
  pull(pos) %>% as.character() %>% unique()
rna_biased_vars_100 <- tf1_data %>% mutate(ratio = (rna_af + pseudo)/(atac_af + pseudo)) %>% 
  filter(ratio > 100 & rna_cov > 100 & atac_cov > 100 & round(rna_af,3) > 0.005 )%>% 
  pull(pos) %>% as.character() %>% unique()
cols <- readRDS("../data/meta_data/cd34_colonyColors.rds")

# Process everything for the two donors
dir <- "../../../mtscATACpaper_large_data_files/source/other/scRNAseq_mgatk_benchmark/"
do_all <- function(donor){
  SE <- readRDS(paste0(dir,"/",donor,"mito.rds"))
  SE <- SE[,colData(SE)$depth > 50]
  
  mut_se <- call_mutations_mgatk(SE)
  misc_df <- data.frame(rowData(mut_se))
  misc_df$old_var <- paste0(as.character(misc_df$position), "_", substrRight1(misc_df$variant))
  saveRDS(mut_se, file = paste0("../../../mtscATACpaper_large_data_files/source/other/",donor,"_smartseq2_SE_af.rds"))
  
  two_df <- read.table("../data/meta_data/SupplementalTable2.tsv", header = TRUE) %>%
    filter(Donor == donor)
  misc_df$known_var <- as.character(misc_df$old_var) %in% as.character(two_df$Var)
  
  # make plot 1 -- new variant calling plot
  p1 <- ggplot(misc_df %>% filter(n_cells_over_5 >= 5 & mean_coverage > 50) %>% 
                 mutate(new_var_color = ifelse(variant %in% c("1672A>G", "2619A>G"), "edited", ifelse(known_var, "known", "unknown"))), 
               aes(x = strand_correlation, y = log10(vmr), color = new_var_color)) +
    geom_point(size = 0.3) + scale_color_manual(values = c("dodgerblue3", "firebrick", "black")) +
    labs(color = "found_prev", x = "Strand concordance", y = "log10 VMR") +
    pretty_plot(fontsize = 8) + L_border() +
    theme(legend.position = "bottom") +
    geom_vline(xintercept = 0.65, linetype = 2) +
    geom_hline(yintercept = -2, linetype = 2) +
    ggtitle(paste0(donor, " SMART-seq2 Colonies"))
  
  cowplot::ggsave2(p1, file = paste0("../plots/", donor, "_varCall.pdf"), width = 2, height = 2.45)
  
  filter_df <- misc_df %>% filter(n_cells_over_5 >= 5) %>% filter(strand_correlation >= 0.65 & log10(vmr) > -2 & mean_coverage > 50)
  
  # Compute QC stats
  old_total = dim(two_df)[1]
  n_recovered <- sum(filter_df %>% pull(known_var))
  rate_hit <- n_recovered/old_total
  n_missed <- old_total - n_recovered
  plot_vars <- filter_df %>% filter(!known_var) %>% pull(variant)
  plot_vars_old_notation <- filter_df %>% filter(!known_var) %>% pull(old_var)
  
  summary_stats_df <- data.frame(donor, n_recovered, n_missed,
                                 new_vars = length(plot_vars), total_new_method = dim(filter_df)[1], rate_hit)
  
  # Look at what we missed
  old_all_df <- misc_df  %>% filter(known_var)
  old_all_df$recovered <- old_all_df$variant %in% as.character(filter_df$variant)
  saveRDS(old_all_df, file = paste0("../output/", donor, "_supervised_variants.rds"))
  saveRDS(filter_df, file = paste0("../output/", donor, "_all_mgatk_variants.rds"))
  
  # Plot the new variants
  afp <- assays(mut_se)[["allele_frequency"]][plot_vars,]
  colony <- str_split_fixed(colnames(afp), "_", 3)[,2]
  
  
  # Make the column annotation for the heatmap
  ha1 <- HeatmapAnnotation(df = data.frame(colony),
                           col = list(colony = cols))
  
  # Thresholds for visualization consistent with last paper
  afp[afp < 0.01] <- 0
  afp[afp > 0.2] <- 0.2
  
  # new found clone-specific variants
  if(donor == "Donor1"){
    clone_vars <- c("1769G>A", "2013G>A", "6506T>C", "2253A>G", "2310A>G", "5110A>G", "11720A>G")
    
  } else {
    clone_vars <- c("2043T>C", "10550T>C", "11344T>C", "1846A>G") 
    
  }
  
  # Add a row annotation for variants that are biased
  rna_vars <- ifelse(plot_vars_old_notation %in% clone_vars, "FALSE", 
                     ifelse(plot_vars_old_notation %in% rna_biased_vars_100, "TRUE2", ifelse(plot_vars_old_notation %in% rna_biased_vars_10, "TRUE1", "FALSE")))
  row_df <- data.frame(
    variant = plot_vars,
    RNA = rna_vars,
    CL = plot_vars %in% clone_vars
  ) %>% arrange(desc(CL), desc(RNA))
  
  ha_row <- rowAnnotation(df = row_df[,c("CL", "RNA"), drop = FALSE], width = unit(0.5, "cm"),
                          col = list(RNA = c("TRUE1"='blue', "TRUE2" = "purple", "FALSE"='white'),
                                     CL = c("TRUE"='blue', "FALSE"='white')))
  
  # visualize new mutations
  pdf(file=paste0("../plots/",donor,"_new_mutations.pdf"), width = 5, height = 3)  
  par(cex.main=0.8,mar=c(1,1,1,1))
  hm <- Heatmap(data.matrix(afp[as.character(row_df$variant),]), col=c("white", as.character(jdb_palette("brewer_red",type="continuous"))),
                cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE,
                row_names_gp = gpar(fontsize = 5),
                top_annotation = ha1,
                name = "AF")
  print(ha_row + hm)
  dev.off()
  summary_stats_df
}

# Make all of the plots necessary for comparing an individual donor
do_all("Donor1")
do_all("Donor2")

# Import previously identified (supervised) variants
d1df <- readRDS("../output/Donor1_supervised_variants.rds"); d1df$donor <- "Donor1"
d2df <- readRDS("../output/Donor2_supervised_variants.rds"); d2df$donor <- "Donor2"

# Set up a wilcox test
wt_df <- rbind(
  d1df, d2df
)
wilcox.test(mean ~ recovered, wt_df)
wt_df$plot_x <- ifelse(wt_df$recovered, "aYes", "No")

# Plot summary of variants replicated or not
p1 <- ggplot(wt_df, aes(x = plot_x, y = mean, color = donor)) +
  geom_quasirandom(size = 0.3) +
  geom_boxplot(color = "dodgerblue3", outlier.shape = NA, fill = NA) +
  scale_color_manual(values = c("black", "grey")) +
  scale_y_log10() + pretty_plot(fontsize = 8) + L_border() +
  labs(x = "Replicated", y = "Population allele frequency (log10)", color = "Donor")

cowplot::ggsave2(p1, file = "../plots/did_not_recover.pdf", width = 2.5, height = 2)
