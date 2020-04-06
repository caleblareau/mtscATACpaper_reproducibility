library(data.table)
library(dplyr)
library(ggrastr)
load("../output/trajectory_inferences.18march2020.rda")

# Import data
dd <- readRDS("../output/cluster-id-800.rds")
dd$timepoint <- as.numeric(factor(str_split_fixed(as.character(dd$barcode), '-', 2)[,1]))
mdf <- merge(dd, df_ery, by.x = "barcode",by.y = "row.names")

# Fine putative clusters of interest
mdf %>% group_by(cluster,Group) %>% summarize(count = n()) %>%
  reshape2::dcast(cluster~Group, value.var = "count", fill = 0) -> odf
odf %>% filter(prog > 2) %>% arrange(desc(my4/(ery4 + ery5 + ery6 + 1)))



# Summarize more values per clone
mdf[complete.cases(mdf),] %>% group_by(cluster, timepoint) %>%
  summarize(count = n(), ps = mean(pseudotime))  %>%
  dplyr::filter(count >= 10) %>% group_by(cluster) %>%
  dplyr::filter(n() > 1) %>% summarize(minps = min(ps),  maxps = max(ps))  %>%
  arrange(desc(maxps)) %>% data.frame() -> mdf2

# Function to extract barcodes from a specific timepoint and a specific cluster
pull_barcodes_dd <- function(timepoint1, cluster1){
  dd %>% dplyr::filter(as.character(cluster) %in% as.character(cluster1)  & as.character(timepoint) %in% as.character(timepoint1) ) %>%
    pull(barcode)
}

plot_d <-  data.frame(
  df_ery, barcode = rownames(df_ery)
)

# Extract interesting 
plot_d$bi1 <- plot_d$barcode  %in% pull_barcodes_dd(1, 30)
plot_d$bi2 <- plot_d$barcode  %in% pull_barcodes_dd(2, 30)
plot_d$bi3 <- plot_d$barcode  %in% pull_barcodes_dd(3,30)
plot_d$ery1 <- plot_d$barcode  %in% pull_barcodes_dd(1, 10)
plot_d$ery2 <- plot_d$barcode  %in% pull_barcodes_dd(2, 10)
plot_d$ery3 <- plot_d$barcode  %in% pull_barcodes_dd(3,10)
plot_d$my1 <- plot_d$barcode  %in% pull_barcodes_dd(1, 54)
plot_d$my2 <- plot_d$barcode  %in% pull_barcodes_dd(2, 54)
plot_d$my3 <- plot_d$barcode  %in% pull_barcodes_dd(3,54)

# find mutations associated with the anecdote clones that were extracted
afin <- assays(readRDS("../output/filtered_mitoSE_CD34-800.rds"))[["allele_frequency"]]

round(Matrix::rowMeans(afin[,as.character(plot_d %>% dplyr::filter(my3) %>% pull(barcode))]) %>% sort(), 3) %>% tail()
round(Matrix::rowMeans(afin[,as.character(plot_d %>% dplyr::filter(ery3) %>% pull(barcode))]) %>% sort(), 3) %>% tail()
round(Matrix::rowMeans(afin[,as.character(plot_d %>% dplyr::filter(bi3) %>% pull(barcode))]) %>% sort(), 3) %>% tail()

# Recode names for convience
plot_d$ery12 <- ifelse(plot_d$ery1, "Ery1", ifelse(plot_d$ery2, 'Ery2', 'A'))
plot_d$my12 <- ifelse(plot_d$my1, "My1", ifelse(plot_d$my2, 'My2', 'A'))
plot_d$bi12 <- ifelse(plot_d$bi1, "Bi1", ifelse(plot_d$bi2, 'Bi2', 'A'))

tb <- theme(legend.position = "none",
            panel.grid = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.background = element_blank()) 


# Visualze clone plots
pe12 <- ggplot(plot_d %>% arrange(ery12), aes(x =  x, y  =  y, color = ery12)) +
  geom_point_rast(raster.dpi = 500, size = 3) + tb + scale_color_manual(values = c("lightgrey", "dodgerblue", "#008080"))

pe3 <- ggplot(plot_d %>% arrange(ery3), aes(x =  x, y  =  y, color = ery3)) +
  geom_point_rast(raster.dpi = 500, size = 3) + tb + scale_color_manual(values = c("lightgrey", "purple3"))

pm12 <- ggplot(plot_d %>% arrange(my12), aes(x =  x, y  =  y, color = my12)) +
  geom_point_rast(raster.dpi = 500, size = 3) + tb + scale_color_manual(values = c("lightgrey", "dodgerblue", "#008080"))

pm3 <- ggplot(plot_d %>% arrange(my3), aes(x =  x, y  =  y, color = my3)) +
  geom_point_rast(raster.dpi = 500, size = 3) + tb + scale_color_manual(values = c("lightgrey", "purple3"))

pb12 <- ggplot(plot_d %>% arrange(bi12), aes(x =  x, y  =  y, color = bi12)) +
  geom_point_rast(raster.dpi = 500, size = 3) + tb + scale_color_manual(values = c("lightgrey", "dodgerblue", "#008080"))

pb3 <- ggplot(plot_d %>% arrange(bi3), aes(x =  x, y  =  y, color = bi3)) +
  geom_point_rast(raster.dpi = 500, size = 3) + tb + scale_color_manual(values = c("lightgrey", "purple3"))

cowplot::ggsave2(cowplot::plot_grid(pe12, pm12, pb12, pe3, pm3, pb3, nrow = 2),
                filename = "../plots/800_raster_clones.png",
                width = 4.2, height = 2.8, units = "in", dpi = 1000)


#------

# Classify the clone type
ery_outcome <- df_ery %>% mutate(barcode = row.names(df_ery)) %>% dplyr::filter(pseudotime >= 50) %>% pull(barcode)
mye_outcome <- df_my %>% mutate(barcode = row.names(df_my)) %>% dplyr::filter(pseudotime >= 50) %>% pull(barcode)
dd$outcome <- ifelse(dd$barcode %in% ery_outcome, "Ery", ifelse(dd$barcode %in% mye_outcome, "Myeloid", "none")) 
dff <- dd %>% dplyr::filter(timepoint == "3") %>%
  group_by(cluster) %>% mutate(count = n()) %>% dplyr::filter(count >= 10)

# Compute the proportion of late timepoint cells observed in each
observed <- dff %>% group_by(cluster) %>%
  summarize(p_ery = sum(outcome == "Ery")/n(), p_myeloid = sum(outcome == "Myeloid")/n())

# Do 100 permutations for a background of assignment proportions
lapply(1:100, function(idx){
  set.seed(idx)
  dff$cluster <- sample(dff$cluster)
  permuted <- dff %>% group_by(cluster) %>%
    summarize(p_ery = sum(outcome == "Ery")/n(), p_myeloid = sum(outcome == "Myeloid")/n(), permutation = idx)
  permuted
}) %>% rbindlist() %>% data.frame() -> permdf


# Compute z-scores using a background permutation of the cluster labels
means_df <- permdf %>% group_by(cluster) %>%
  summarize(m_ery = mean(p_ery), sd_ery = sqrt(var(p_ery)), 
            m_mye = mean(p_myeloid), sd_mye = sqrt(var(p_myeloid)))

observed$eryZ <- (observed$p_ery - means_df$m_ery)/means_df$sd_ery
observed$myeZ <- (observed$p_myeloid - means_df$m_mye)/means_df$sd_mye
observed$diff <- observed$eryZ - observed$myeZ
observed <- observed %>% arrange(desc(diff)) %>% mutate(rankSort = 1:n())


# Plot rank-sorted z-scores; cap at 20 for visualization purposes
p1 <- ggplot(observed, aes(x = rankSort, y = ifelse(diff >20, 20, diff), fill = ifelse(diff >20, 20, diff))) +
  geom_bar(stat = "identity", color = "black") +
  pretty_plot(fontsize = 8) + L_border() + 
  geom_hline(yintercept = 0) +
  theme(legend.position = "none") +
  scale_fill_gradientn(colors = jdb_palette("brewer_spectra")) +
  labs(x = "Rank sorted clones", y = "Erythroid Z - Myeloid Z")

cowplot::ggsave2(p1, file = "../plots/zscore_clones_bars.pdf", height = 1.5, width = 2.5)

# Infer the biased clones
ery_biased_clones <- observed %>% dplyr::filter(diff >= 5) %>% pull(cluster) %>% as.character()
mye_biased_clones <- observed %>% dplyr::filter(diff <= -5) %>% pull(cluster) %>% as.character()

length(mye_biased_clones)
length(ery_biased_clones)

# pull the barcodes from the first timepoint that were associated with the biased clones
ery_bc <- dd %>% dplyr::filter(timepoint == "1" & cluster %in% ery_biased_clones) %>% pull(barcode) %>% as.character()
mye_bc <- dd %>% dplyr::filter(timepoint == "1" & cluster %in% mye_biased_clones) %>% pull(barcode) %>% as.character()

# Subset to the progenitors for followup analysis
if(FALSE){
  prog <- df_ery %>% mutate(barcode = rownames(df_ery)) %>% dplyr::filter(Group == "prog") %>% pull(barcode)
  ery_bc <- intersect(prog, ery_bc)
  mye_bc <- intersect(prog, mye_bc)
}

# Import deviation scores from chromVAR to see what TFs are 
dev <- assays(readRDS("../output/chromVAR/CD34_800_Day08_tf_deviations.rds"))[["z"]]
colnames(dev) <- paste0("CD34_800_Day08-", as.character(colnames(dev)))
dev_df <- data.frame(
  rowData(readRDS("../output/chromVAR/CD34_800_Day08_tf_deviations.rds")),
  ery_biased = rowMeans(dev[,ery_bc]),
  mye_biased = rowMeans(dev[,mye_bc])
)
dev_df2 <- dev_df %>% mutate(diff = ery_biased - mye_biased) %>%
  arrange(desc(diff)) %>% mutate(rank = 1:n())

# Export pretty plot
p2 <- ggplot(dev_df2, aes(x = rank, y = diff, color = diff)) +
  scale_color_gradientn(colors = jdb_palette("brewer_spectra")) +
  pretty_plot(fontsize = 8) + L_border() + 
  geom_hline(yintercept = 0) +
  theme(legend.position = "none") +
  geom_point(size = 0.5) +
  labs(x = "Rank sorted TF motifs", y = "Difference in deviations (z-score)")
cowplot::ggsave2(p2, file = "../plots/zscore_TFs_dots.pdf", height = 1.5, width = 2.5)


# Visualize plot that LL had in mind
prog <- df_ery %>% mutate(barcode = rownames(df_ery)) %>% dplyr::filter(Group == "prog")
prog$color <- ifelse(prog$barcode %in% ery_bc, "erythroid",
                     ifelse(prog$barcode %in% mye_bc, "myeloid",
                            ifelse(substr(prog$barcode, 6,6) == "5", "500culture", "balanced")))

ggplot(prog %>% filter(color %in% c("erythroid", "myeloid")), aes(x = x, y = y, color = color)) +
  geom_point() + scale_color_manual(values = c("firebrick", "orange"))
