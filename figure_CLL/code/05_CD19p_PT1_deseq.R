library(BuenColors)
library(stringr)
library(SummarizedExperiment)
library(dplyr)
library(data.table)
library(Matrix)
library(DESeq2)
library(ggrastr)
register(MulticoreParam(2))
"%ni%" <- Negate("%in%")
source("00b_mgatk_helper_fn.R")

SE <- readRDS("../../../mtscATACpaper_large_data_files/source/mgatk_output/CLL1_CD19pos_5p_scRNAseq_mgatk.rds")

# Assemble data frame of useful mutations
odf <- data.frame(
  Barcode = colnames(SE),
  e2(5140, "G", "A"),
  e2(1260, "A", "G"),
  e2(14858, "G", "A"),
  e2(1872, "T", "C")
)

# Look at the number of cells per 'clone'
odf %>% filter(het5140G_A > 0.9 & cov5140 >= 1) %>% dim()
odf %>% filter(het1260A_G > 0.9 & cov1260 >= 1) %>% dim()
odf %>% filter(het1872T_C > 0.9 & cov1872 >= 1) %>% dim()
odf %>% filter(het14858G_A > 0.9 & cov14858 >=1) %>% dim()

# Pull barcodes of near-homoplasmic mutations
bc5140 <- odf %>% filter(het5140G_A > 0.9) %>% pull("Barcode") %>% as.character()
bc1872 <- odf %>% filter(het1872T_C > 0.9) %>% pull("Barcode") %>% as.character()
bc14858 <- odf %>% filter(het14858G_A > 0.9) %>% pull("Barcode") %>% as.character()
bc1260 <- odf %>% filter(het1260A_G > 0.9) %>% pull("Barcode") %>% as.character()

# Import genes
dt <- fread("zcat < ../scRNAseq_data/gene_expression/CLL1_CD19pos_ge/matrix.mtx.gz", skip = 3)
barcodes <- fread("../scRNAseq_data/gene_expression/CLL1_CD19pos_ge/barcodes.tsv", header = FALSE)[[1]]
genes <- make.unique(fread("../scRNAseq_data/gene_expression/CLL1_CD19pos_ge/features.tsv.gz", header = FALSE)[[2]])

# Assemble matrix
mat_pos <- sparseMatrix(
  i = c(dt[["V1"]],length(genes)),
  j = c(dt[["V2"]],1),
  x = c(dt[["V3"]], 0)
)
rownames(mat_pos) <- genes; colnames(mat_pos) <- barcodes

# Compute all genes that are differential between two sets of barcodes
compute_pairwise <- function(bc1, bc2){
  
  # Subset to relevant genes
  mat_ss <- mat_pos[,c(bc1, bc2)]
  rs <- rowSums(mat_ss)
  boo <- (rs/sum(rs) * 1000000) > 1 & !grepl("^MT-", rownames(mat_ss)) # filter mito genes and genes with a count less than 1 cpm
  RNA.counts <- data.matrix(mat_ss[boo,])
  
  #Setup DEseq2
  RNA.counts.df <- as.data.frame(RNA.counts)
  
  # Establish column data
  RNA.condition <- c(rep("bc1", length(bc1)), rep("bc2", length(bc2)))
  colData <- as.data.frame(RNA.condition)
  row.names(colData) <- colnames(RNA.counts.df)
  
  # Run DEseq2
  RNA.dds <- DESeqDataSetFromMatrix(countData = RNA.counts.df, colData = colData, design = ~ RNA.condition)
  RNA.dds <- DESeq(RNA.dds, parallel = TRUE)
  resSig.RNA <- data.frame(results(RNA.dds, contrast=c("RNA.condition", "bc1", "bc2"), 
                                   parallel = TRUE, alpha = 0.01, independentFiltering = FALSE))
  resSig.RNA$gene <- rownames(RNA.counts)
  resSig.RNA <- resSig.RNA %>% arrange(pvalue)
  resSig.RNA
}

# Just to observe the number of differentially expressed genes
compute_pairwise(bc14858, bc5140) %>% filter(padj < 0.05) %>% dim() # 11
compute_pairwise(bc14858, bc1872) %>% filter(padj < 0.05) %>% dim() # 12
compute_pairwise(bc14858, bc1260) %>% filter(padj < 0.05) %>% dim() # 6

compute_pairwise(bc1872, bc5140) %>% filter(padj < 0.05) %>% dim() # 4
compute_pairwise(bc1872, bc1260) %>% filter(padj < 0.05) %>% dim() # 1


# Do comparison between the major clone and the others
resSig.RNA <- compute_pairwise(bc14858, unique(c(bc1260, bc5140, bc1872))) 
resSig.RNA$color <- ifelse(resSig.RNA$padj < 0.05, ifelse(resSig.RNA$gene %in% c("IGLV3-1", "IGKV3-11", "IGLV2-5", "IGLC3", "IGLC2", "IGHV1-69D"), "heavy", "anno"), "other")

# Cap the results for plotting purposes
resSig.RNA$log2FoldChange <- ifelse(resSig.RNA$log2FoldChange > 4, 4, resSig.RNA$log2FoldChange)
resSig.RNA$log2FoldChange <- ifelse(resSig.RNA$log2FoldChange < -4, -4, resSig.RNA$log2FoldChange)
resSig.RNA$padj <- pmax(1e-10, resSig.RNA$padj)

# Visualize in a volcano
p1 <- ggplot(resSig.RNA, aes(x = log2FoldChange, y = -log10(padj), color = color)) +
  geom_point_rast(size = 0.4, raster.dpi = 1000) + 
  geom_text(data=subset(resSig.RNA, padj < 0.05 & color != "heavy"), size = 1.5,
            aes(log2FoldChange,-log10(padj),label=gene), vjust = -0.5) +
  geom_hline(yintercept = 1.3, linetype = 2) +
  pretty_plot(fontsize = 7) + L_border() + scale_color_manual(values = c("dodgerblue3", "purple3" , "black")) +
  scale_x_continuous(limits = c(-4, 4)) + theme(legend.position = "none") +
  labs(x = "log2 minor/major")

cowplot::ggsave2(p1, file = "../plots/CLL_volcano_pt1.pdf", width = 1.7, height = 1.7)


