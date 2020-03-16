library(BuenColors)
library(stringr)
library(SummarizedExperiment)
library(dplyr)
library(data.table)
library(Matrix)
library(ggrastr)
source("00a_signature_function.R")
source("00b_mgatk_helper_fn.R")

# Import mutations
SE <- readRDS("../../../mtscATACpaper_large_data_files/source/mgatk_output/CLL1_CD19neg_5p_scRNAseq_mgatk.rds")

# Pull relevant mutations
odf <- data.frame(
  barcode = colnames(SE),
  e2(5140, "G", "A"),
  e2(1260, "A", "G"),
  e2(14858, "G", "A"),
  e2(1872, "T", "C")
)

# Annotate cell state based on mutation hierarchy
odf$mutation_id <- ifelse(odf$het5140G_A > 0.9 & odf$cov5140 >= 5, "5140G>A", ifelse(odf$het1260A_G > 0.9 & odf$cov1260 >= 5, "1260A>G",
                                                                  ifelse(odf$het14858G_A > 0.9 & odf$cov14858 >= 5, "14858G>A", ifelse(odf$het1872T_C > 0.9 & odf$cov1872 >= 5, "1872T>C",""))))

# Import gene expression counts
dt <- fread("zcat < ../scRNAseq_data/gene_expression/CLL1_CD19neg_ge/matrix.mtx.gz", skip = 3)
barcodes <- fread("../scRNAseq_data/gene_expression/CLL1_CD19neg_ge/barcodes.tsv", header = FALSE)[[1]]
genes <- make.unique(fread("../scRNAseq_data/gene_expression/CLL1_CD19neg_ge/features.tsv.gz", header = FALSE)[[2]])
mat_neg <- sparseMatrix(
  i = c(dt[["V1"]],length(genes)),
  j = c(dt[["V2"]],1),
  x = c(dt[["V3"]], 0)
)
rownames(mat_neg) <- genes; colnames(mat_neg) <- barcodes

# Generate single-cell signatures
module_scores <- project_mats_do_all(mat_neg)
merge_df <- merge(module_scores, odf, by = "barcode")

# Visualize module scores
p1a <- ggplot(merge_df %>% arrange(mutation_id), aes(x = Bcell_signature, y =Tcell_signature, color = mutation_id)) +
  geom_point_rast(size = 0.2, raster.dpi = 1000) + pretty_plot(fontsize = 7) + L_border() +
  scale_color_manual(values = c("lightgrey", "dodgerblue4", "dodgerblue", "orange2","firebrick")) + 
  theme(legend.position = "none") + labs(x = "Bcell Signature", y = "Tcell Signature")

cowplot::ggsave2(p1a, file = "../plots/signatures_CD19neg_CLL1_compressed.pdf", 
                width = 1.7, height = 1.7)


