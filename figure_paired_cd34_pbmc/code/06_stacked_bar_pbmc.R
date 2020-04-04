library(data.table)
library(dplyr)
library(SummarizedExperiment)
library(BuenColors)

sdf <- readRDS('../output/PBMC_clone_DF.rds')

# Make proportions of each lineage by clone
for_p <- sdf %>%
  dplyr::filter(lineage_assignment_atac != "other") %>% 
  group_by(mito_cluster, lineage_assignment_atac) %>% summarize(count = n()) %>%
  ungroup() %>% group_by(mito_cluster) %>% mutate(total = sum(count)) %>%
  mutate(proportion = count/total) %>%
  dplyr::filter(total >= 10) %>% arrange(desc(proportion))

for_p$ca <- factor(as.character(for_p$mito_cluster), levels = unique(as.character(for_p$mito_cluster)))
length(unique(for_p$ca ))

p1 <- ggplot(for_p, aes(fill=lineage_assignment_atac, y=proportion, x=ca)) + 
  pretty_plot(fontsize = 8) +
  geom_bar(position="stack", stat="identity") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_fill_manual(values = c("purple3", "orange2", "dodgerblue3")) +
   L_border() + labs(fill = "") +
  theme(legend.position = "bottom") +
  scale_y_continuous(limits = c(0,1), expand=c(0,0))

cowplot::ggsave2(p1, file = "../plots/filled_bar_PBMCs.pdf", width = 1.8, height = 1.8)

