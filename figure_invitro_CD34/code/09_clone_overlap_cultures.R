library(stringr)
library(data.table)
library(dplyr)
library(matchingR)

load("../output/trajectory_inferences.18march2020.rda")

# Import all variants
se500 <- readRDS("../output/filtered_mitoSE_CD34-500.rds")
se800 <- readRDS("../output/filtered_mitoSE_CD34-800.rds")
clusters_500 <- readRDS("../output/cluster-id-500.rds")
clusters_800 <- readRDS("../output/cluster-id-800.rds")


# Find the variants that are overlapping in the cultures-- subset only to clusters that utilize those mutations at a mean 5%
shared_variants <- intersect(rownames(se800), rownames(se500))

# Extract shared mutations and the clusters that participate in them for each of the cultures
df_500 <- data.frame(
  cluster = clusters_500[,"cluster"], t(data.matrix(assays(se500)[["allele_frequency"]][shared_variants,]))
) %>% group_by(cluster) %>% summarise_all(mean) %>% mutate_if(is.numeric, round, 3)
boo500 <- rowSums(data.matrix(df_500[,-1]) >= 0.05) > 0
cl500mat <- t(data.matrix(df_500[boo500,-1]))
colnames(cl500mat) <- paste0("cl500_", as.character(df_500[[1]]))[boo500]

# Same thing for 800
df_800 <- data.frame(
  cluster = clusters_800[,"cluster"], t(data.matrix(assays(se800)[["allele_frequency"]][shared_variants,]))
) %>% group_by(cluster) %>% summarise_all(mean) %>% mutate_if(is.numeric, round, 3)
boo800 <- rowSums(data.matrix(df_800[,-1]) >= 0.05) > 0
cl800mat <- t(data.matrix(df_800[boo800,-1]))
colnames(cl800mat) <- paste0("cl800_", as.character(df_800[[1]]))[boo800]

# Run optimal matching
cormat <- cor(cl800mat,cl500mat, method = "pearson")
clone_matches <- galeShapley.collegeAdmissions(studentUtils = cormat, collegeUtils = t(cormat), slots = 1)

optim_df <- data.frame(
  clone500 = colnames(cl500mat)[clone_matches$matched.colleges],
  clone800 = colnames(cl800mat)
)

# Now compute the initial clone size for each culture
boo_d8_c5 <- substr(as.character(clusters_500$barcode), 14,14) == 8
prop_5 <- clusters_500[boo_d8_c5,] %>%
  group_by(cluster) %>% summarize(prop = n()/sum(boo_d8_c5))

boo_d8_c8 <- substr(as.character(clusters_800$barcode), 14,14) == 8
prop_8 <- clusters_800[boo_d8_c8,] %>%
  group_by(cluster) %>% summarize(prop = n()/sum(boo_d8_c8))

# Make named vectors to append to optim df
p5v <- prop_5$prop; names(p5v) <- paste0("cl500_", as.character(prop_5[[1]]))
p8v <- prop_8$prop; names(p8v) <- paste0("cl800_", as.character(prop_8[[1]]))

optim_df$prop5 <- p5v[as.character(optim_df[,"clone500"])]
optim_df$prop8 <- p8v[as.character(optim_df[,"clone800"])]

# Now plot the result
pCo <- ggplot(optim_df, aes(x = prop5*100, y = prop8*100)) +
  geom_point(size = 0.2) +
  geom_smooth(method = "lm", se = FALSE, size = 0.5,color = "purple") +
  geom_abline(intercept = 0, slope = 1, linetype = 2, size = 0.5) +
  pretty_plot(fontsize = 8) + L_border() +
  labs(x = "% cells at day 8 (500 input)", y = "% cells at day 8 (800 input)") 

summary(lm(prop8 ~ prop5, optim_df))

cowplot::ggsave2(pCo, file = "../plots/plot_clone_overlap_size.pdf", width = 1.7, height = 1.7)