library(BuenColors)
library(stringr)
library(dplyr)
library(data.table)

# import cells called by cellranger
cells1 <- read.table("../scRNAseq_data/gene_expression/CLL1_CD19pos_ge/barcodes.tsv", header = FALSE, sep = "\t")[,1] %>% as.character()
cells2 <- read.table("../scRNAseq_data/gene_expression/CLL2_CD19pos_ge/barcodes.tsv", header = FALSE, sep = "\t")[,1] %>% as.character()

# Import clonotypes
clonotypes1 <- fread("../scRNAseq_data/BCR_clonotypes/CLL_PT1_BCR_filtered_contig_annotations.csv") %>%
  data.frame() %>% filter(raw_consensus_id != "None") %>% filter(barcode %in% cells1) %>%
  group_by(barcode) %>% top_n(1, wt = "umis")

clonotypes2 <- fread("../scRNAseq_data/BCR_clonotypes/CLL_PT2_BCR_filtered_contig_annotations.csv") %>%
  data.frame() %>% filter(raw_consensus_id != "None") %>% filter(barcode %in% cells2) %>%
  group_by(barcode) %>% top_n(1, wt = "umis")

# Compute rates of most abundant clonotype
pt1_rate <- sum(clonotypes1$raw_clonotype_id == "clonotype1")/length(clonotypes1$raw_clonotype_id) *100
pt2_rate <- sum(clonotypes2$raw_clonotype_id == "clonotype1")/length(clonotypes2$raw_clonotype_id) *100

# Plot rates of most abundant clonotype
p1 <- data.frame(rates = c(pt2_rate, pt1_rate), patient = c("PatientX2", "PatientY1")) %>%
  ggplot(aes(x = patient, y = rates)) +
  geom_bar(stat = "identity", color = "black", fill = "lightgrey", width = 0.7 ) + coord_flip() +
  scale_y_continuous(limits = c(0,100), expand = c(0,0)) +
  pretty_plot(fontsize = 7) + L_border() +
  labs(x = "", y = "% CD19+ cells with major clonotype")
cowplot::ggsave2(p1, file = "../plots/BCRrate_plot.pdf", width = 2, height = 1.7)
