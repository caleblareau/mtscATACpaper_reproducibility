library(BuenColors)
library(ggrastr)
library(SummarizedExperiment)
library(dplyr)

# Import data
clone_df_in <- readRDS("../output/CD34_clone_DF.rds")

# Look at chi-square associations between clones and the cell type cluster
# This is for the CD34 data
get_chisquare_stats <- function(clone_df_input, seed = 0){
  
  if(seed == 0){
    clone_df <- clone_df_input
  } else {
    set.seed(seed)
    clone_df_input$mito_cluster <- sample(clone_df_input$mito_cluster)
    clone_df <- clone_df_input
  }
  
  # Pull the clone abundance
  clones <- clone_df %>% group_by(mito_cluster) %>% summarize(n=n()) %>% dplyr::filter(n >= 10) %>% pull(mito_cluster)
  chromatin_clusters <- as.character(sort(unique(clone_df$Clusters)))
  n_global <- clone_df %>% group_by(Clusters) %>% summarize(prop = n() / dim(clone_df)[1], n = n()) %>% pull(n)
  prop_global <- n_global/sum(n_global)
  
  # Loop over the clones and compute statistics
  lapply(clones, function(clone){
    
    ndf <- clone_df %>% filter(mito_cluster == clone) %>%
      group_by(Clusters) %>% summarize(n = n())
    
    # Make a named vector of observations per cluster
    n_vec <- ndf$n; names(n_vec) <- ndf$Clusters
    n_clone <- unname(n_vec[chromatin_clusters])
    n_clone <- ifelse(is.na(n_clone), 0, n_clone)
    
    cs <- chisq.test(n_clone, p = prop_global)
    
    data.frame(
      clone_id = clone,
      CS_stat = unname(cs$statistic),
      CS_pvalue = unname(cs$p.value),
      n_cells = sum(ndf$n)
    )
  }) %>% rbindlist() %>% data.frame() -> clone_stat_df
  clone_stat_df$FDR <- p.adjust(clone_stat_df$CS_pvalue, "fdr")
  clone_stat_df$log10FDR <- -1*log10(clone_stat_df$FDR)
  clone_stat_df$log10pvalue <- -1*log10(clone_stat_df$CS_pvalue)
  
  print(length(clones))
  clone_stat_df %>% arrange(desc(log10pvalue)) %>% mutate(rank = 1:n(), seed = seed)
}

# Do it for both observed and permuted
cs_obs <- get_chisquare_stats(clone_df_in); cs_obs$what <- "Observed"
cs_perm <- get_chisquare_stats(clone_df_in, 1); cs_perm$what <- "Permuted"

p1 <- ggplot(rbind(cs_obs, cs_perm) %>% arrange(desc(what)), aes(x = rank, y = log10FDR, color = what)) + 
  geom_point(size = 0.4)+
  scale_color_manual(values = c("black", "lightgrey"))  + 
  labs(x = "Rank sorted clones", y = "-log10 FDR") + 
  pretty_plot(fontsize = 7) + L_border() +  scale_y_continuous(limits = c(0,2)) +
  theme(legend.position = "none") + ggtitle("CD34+ mtDNA clones")


# Harmonize the pbmc data a bit
pbmc_clone_df_in <- readRDS("../output/PBMC_clone_DF.rds")
pbmc_clone_df_in <- pbmc_clone_df_in %>% dplyr::filter(lineage_assignment_atac != "other")
pbmc_clone_df_in$Clusters <- pbmc_clone_df_in$lineage_assignment_atac

cs_obs_pbmc <- get_chisquare_stats(pbmc_clone_df_in); cs_obs_pbmc$what <- "Observed"
cs_perm_pbmc <- get_chisquare_stats(pbmc_clone_df_in, 1); cs_perm_pbmc$what <- "Permuted"

p2 <- ggplot(rbind(cs_obs_pbmc, cs_perm_pbmc) %>% arrange(desc(what)), aes(x = rank, y = log10FDR, color = what)) + 
  geom_point(size = 0.4)+
  scale_color_manual(values = c("black", "lightgrey"))  + 
  labs(x = "Rank sorted clones", y = "-log10 FDR") + 
  pretty_plot(fontsize = 7) + L_border() + scale_y_continuous(limits = c(0,2)) +
  theme(legend.position = "none") + ggtitle("PBMC mtDNA clones")

cowplot::ggsave2(cowplot::plot_grid(p1, p2, nrow = 1), height = 1.7, width = 3.4, file = "../plots/bias_X2_comparison.pdf")

head(cs_obs)
head(cs_obs_pbmc)

# Look at the fold change of things
source_df <- data.frame(cluster = c(clone_df_in$mito_cluster, pbmc_clone_df_in$mito_cluster),
                        source = c(rep("CD34", length(clone_df_in$mito_cluster)),
                                   rep("PBMC", length(pbmc_clone_df_in$mito_cluster)))) %>%
  group_by(cluster,source) %>% summarize(count = n())

dd <- reshape2::dcast(source_df, cluster ~ source, value.var = "count", fill = 1)
dd$CD34p <- dd$CD34/sum(dd$CD34)
dd$PBMCp <- dd$PBMC/sum(dd$PBMC)

observed_foldchange <- log2((dd$PBMCp )/(dd$CD34p))

set.seed(1)
perm_df <- data.frame(cluster = c(clone_df_in$mito_cluster, pbmc_clone_df_in$mito_cluster),
                      source = sample(c(rep("CD34", length(clone_df_in$mito_cluster)),
                                        rep("PBMC", length(pbmc_clone_df_in$mito_cluster))))) %>%
  group_by(cluster,source) %>% summarize(count = n())

dd2 <- reshape2::dcast(perm_df, cluster ~ source, value.var = "count", fill = 1)
dd2$CD34p <- dd2$CD34/sum(dd2$CD34)
dd2$PBMCp <- dd2$PBMC/sum(dd2$PBMC)

permuted_foldchange <- log2((dd2$PBMCp )/(dd2$CD34p))

density_df <- data.frame(
  logFC = c(observed_foldchange, permuted_foldchange),
  what = c(rep("observed", length(observed_foldchange)),
               rep("permuted", length(permuted_foldchange))
  ))
pDens <- ggplot(density_df, aes(x = logFC, fill = what)) + 
  geom_density(alpha = 0.2) +
  scale_fill_manual(values =c ("firebrick", "grey")) +
  scale_y_continuous(expand = c(0,0)) +
  pretty_plot(fontsize = 8) + L_border() +
  labs(x = "log2 FC", y = "empirical density") +
  theme(legend.position = "none")
cowplot::ggsave2(pDens, file = "../plots/CD34_PBMC_density.pdf", width = 1.7, height = 1.7)

ks.test(permuted_foldchange, observed_foldchange)
