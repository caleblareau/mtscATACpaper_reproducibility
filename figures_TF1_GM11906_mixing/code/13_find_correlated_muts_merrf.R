library(SummarizedExperiment)
library(Matrix)
library(dplyr)
library(BuenColors)
library(ComplexHeatmap)
library(circlize)
library(ggrastr)

"%ni%" <- Negate("%in%")
perc.rank <- function(x) trunc(rank(x))/length(x)

SE <- cbind(readRDS("../data/big_gi/Mix_Fix_1h_CR-mtMask_mgatk.rds"), 
            readRDS("../data/big_gi/Mix_Fix_6h_CR-mtMask_mgatk.rds"))
SE_chrom <- readRDS("../data/big_gi/GM11906_combinedSE.rds")
df <- data.frame(colData(SE_chrom))

SE <- SE[,as.character(df$cell_id)]
SE_0 <- SE[,as.character(df$cell_id)[df$merrf_af == 0]]

# Process mutations
cov <- assays(SE)[["coverage"]]+ 0.001
ref_allele <- as.character(rowRanges(SE)$refAllele)

getMutMatrix <- function(letter){
  mat <- (assays(SE)[[paste0(letter, "_counts_fw")]] + assays(SE)[[paste0(letter, "_counts_rev")]]) / cov
  rownames(mat) <- paste0(as.character(1:dim(mat)[1]), toupper(ref_allele), ">", letter)
  return(mat[toupper(ref_allele) != letter,])
}

rowVars <- function(x, ...) {
  rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}

mmat <- rbind(getMutMatrix("A"), getMutMatrix("C"), getMutMatrix("G"), getMutMatrix("T"))

plot_df <- data.frame(m8344 = mmat["8344A>G",], m8202 = mmat["8202T>C",], colData(SE)) 

if(FALSE){
  write.table(plot_df[plot_df$m8202 > 0.05 & plot_df$m8344 > 0.05,], file = "../output/double_positive_population.tsv", 
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}
plot_df$color <- perc.rank(plot_df$depth)

# Make plot of scatter between two mutations
p1 <- ggplot(plot_df, aes(x = m8344*100, y = m8202*100, color = color)) + 
  geom_point_rast(dpi = 2000, size = 0.4) +
  scale_color_gradientn(colors = jdb_palette("solar_basic")) +
  pretty_plot(fontsize = 8) + L_border() + labs (x = "m.8344A>G (%)", y = "m.8202T>C (%)") +
  theme(legend.position = "none")

cowplot::ggsave(p1, file = "../output/plots/co-8202_8344.pdf", width = 1.7, height = 1.7)


df <-  data.frame(
  cor = cor(data.matrix(mmat["8344A>G",]), t(data.matrix(mmat)), method = "spearman")[1,],
  mutation = rownames(mmat)
) 

# Show all mutations in a correlation plot

cord_plot_df <- df[complete.cases(df), ] %>% arrange(desc(cor)) %>%
  filter(mutation != "8344A>G") %>% 
  mutate(rank = 1:n())

pA <- ggplot(cord_plot_df, aes(x = rank, y = cor)) +
  geom_point_rast(dpi = 2000, size = 0.2) + pretty_plot(fontsize = 8) +
  L_border() + labs(x = "Rank", y = "Correlation with m.8344A>G") +
  geom_hline(yintercept = 0, linetype = 2) + scale_x_log10()
cowplot::ggsave(pA, file = "../output/plots/correlation_8344.pdf", width = 1.7, height = 1.7)
