library(BuenColors)
library(stringr)
library(SummarizedExperiment)
library(dplyr)
library(data.table)
library(Matrix)
library(ComplexHeatmap)
library(circlize)
register(MulticoreParam(2))
"%ni%" <- Negate("%in%")
source("00b_mgatk_helper_fn.R")

SE <- readRDS("../../../mtscATACpaper_large_data_files/source/mgatk_output/CLL1_CD19pos_5p_scRNAseq_mgatk.rds")


# Produce a per-barcode summary of mutations
odf <- data.frame(
  barcode = colnames(SE),
  e2(5140, "G", "A"),
  e2(1260, "A", "G"),
  e2(14858, "G", "A"),
  e2(1872, "T", "C")
)


# Import and filter for relevant cells the BCR clonotypes
cells1 <- read.table("../scRNAseq_data/gene_expression/CLL1_CD19pos_ge/barcodes.tsv", header = FALSE, sep = "\t")[,1] %>% as.character()
clonotypes1 <- fread("../scRNAseq_data/BCR_clonotypes/CLL_PT1_BCR_filtered_contig_annotations.csv") %>%
  data.frame() %>% filter(raw_consensus_id != "None") %>% filter(barcode %in% cells1) %>%
  group_by(barcode) %>% top_n(1, wt = "umis")

#Filter for cells that have at least a coverage of 2 for one of the important mutation
merge_df <- merge(odf, clonotypes1, by.x = "barcode", by.y = "barcode", all.x = FALSE, all.y = FALSE)
short_df <- merge_df[(merge_df$cov5140 >= 2 | merge_df$cov1260 >= 2 | merge_df$cov14858 >= 2 | merge_df$cov1872 >= 2),c("het14858G_A", "het1260A_G", "het5140G_A", "het1872T_C", "raw_clonotype_id")] 

# Remove doublets from both clonotypes and mtDNA-- clonotypes 1 and 5 are likely the true 'major' clonotypes (see html report); clones 2-4 are likely some sort of cell doublet
rs <- rowSums(data.matrix(short_df[,c(1:4)]))
out_df <- short_df[((rs > 0.99 & rs < 1.01)) & short_df$raw_clonotype_id %ni% c("clonotype2", "clonotype3"),]
out_df$ctid <- ifelse(out_df$raw_clonotype_id %in% c("clonotype1", "clonotype5"), out_df$raw_clonotype_id , "other")
out_df$score <- round(out_df$het14858G_A)*1 +  round(out_df$het1260A_G)*2 + round(out_df$het5140G_A)*3 + round(out_df$het1872T_C)*4
out_df <- out_df %>% arrange(desc(score), (ctid))

# Assign colors to clonotypes
color_vec <- c("dodgerblue2", "purple3", "black"); names(color_vec) <- c("clonotype1", "clonotype5", "other")
ha_col <- HeatmapAnnotation(df = data.frame(Clonotype = out_df$ctid),
                            col = list(Clonotype = color_vec))

# Visualize heatmap
pdf(paste0("../plots/PT1_scRNA_BCR_mito_clonotype.pdf"), width=4, height=2)
hm <- Heatmap(t(data.matrix(out_df[,c(4:1)])), 
              col=as.character(jdb_palette("solar_rojos",type="continuous")),
              show_row_names = TRUE, 
              cluster_columns = FALSE,
              top_annotation=ha_col,
              name = "AF",
              cluster_rows = FALSE, 
              show_column_names = FALSE)
hm
dev.off()


