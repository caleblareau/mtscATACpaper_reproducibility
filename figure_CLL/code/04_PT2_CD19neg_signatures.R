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
SE <- readRDS("../../../mtscATACpaper_large_data_files/source/mgatk_output/CLL2_CD19neg_5p_scRNAseq_mgatk.rds")

odf <- data.frame(
  barcode = colnames(SE),
  e2(4853, "G", "A"),
  e2(12980, "G", "A")
)

# Iport the nuclear mutations from hisnpper
nm <- fread("../scRNAseq_data/nuclear_mutations_pt2/hisnipper_CD19neg.tsv")
odf$LEF1 <- odf$barcode %in% (nm %>% filter(V3 == "C" & V2 == 109084804) %>% pull(V5) %>% unique())
odf$HCST <- odf$barcode %in% (nm %>% filter(V3 == "A" & V2 == 36394730) %>% pull(V5) %>% unique())

# be slightly more strict about coverage (10x) because there's a ton of mito mutations otherwise, mostly in B-cells, and it gets over plotted
odf$score <- ifelse(odf$het4853G_A > 0.9 & odf$cov4853 >= 10, "4853G>A", ifelse(odf$het12980G_A > 0.9 & odf$cov12980 >= 10, "12980G>A", 
                                                                                        ifelse(odf$HCST, "HCST", ifelse(odf$LEF1, "LEF1",ifelse(odf$het4853G_A > 0.9, "4853G>A", ifelse(odf$het12980G_A > 0.9, "12980G>A","znone"))))))

# Import gene expression counts
dt <- fread("zcat < ../scRNAseq_data/gene_expression/CLL2_CD19neg_ge/matrix.mtx.gz", skip = 3)
barcodes <- fread("../scRNAseq_data/gene_expression/CLL2_CD19neg_ge/barcodes.tsv", header = FALSE)[[1]]
genes <- make.unique(fread("../scRNAseq_data/gene_expression/CLL2_CD19neg_ge/features.tsv.gz", header = FALSE)[[2]])
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
p1a <- ggplot(merge_df %>% arrange(desc(score)), aes(x = Bcell_signature, y =Tcell_signature, color = score)) +
  geom_point_rast(size = 0.2, raster.dpi = 1000) + pretty_plot(fontsize = 7) + L_border() +
  scale_color_manual(values = c("znone"="lightgrey", "LEF1" = "dodgerblue4","HCST" = "orange2", "4853G>A" = "dodgerblue", "12980G>A"="firebrick")) + 
  theme(legend.position = "none") + labs(x = "Bcell Signature", y = "Tcell Signature")

cowplot::ggsave2(p1a, file = "../plots/signatures_CD19neg_CLL2_compressed.pdf", 
                width = 1.7, height = 1.7)

