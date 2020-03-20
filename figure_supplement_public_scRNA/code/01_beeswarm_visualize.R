library(data.table)
library(dplyr)
library(ggbeeswarm)

importTop <- function(file, what){
  nt = 500
  cov <- fread(file, header = FALSE) %>% data.frame() %>% 
    arrange(desc(V2)) %>% head(nt) %>% pull(V2)
  df <- data.frame(
    Technology = what, 
    Coverage = cov
  )
  df
}

df <- rbind(
  importTop("../data/depthtables/pbmc-v2-8k-mgatk.depthTable.txt", "10x_V2"),
  importTop("../data/depthtables/pbmc-v3p-10k-mgatk.depthTable.txt", "10x_V3"),
  importTop("../data/depthtables/pbmc-v5p-8k-mgatk.depthTable.txt", "10x_5p"),
  importTop("../data/depthtables/smartseq2-colonies-mgatk.depthTable.txt", "SmartSeq2")
)

p1 <- ggplot(df, aes(x = Technology, y = Coverage)) +
  geom_quasirandom_rast(raster.dpi = 1000, size = 0.5) +
  pretty_plot(fontsize = 8) + L_border() +
  labs(x = "Technology (Top 500 cells / Library)", y = "Coverage") +
  theme(legend.position = "bottom")  +
  scale_y_log10(breaks = c(10, 50, 100, 500, 1000), limits = c(9, 2000))

cowplot::ggsave2(p1, file = "../plots/comparison_technologies_paper.pdf", width = 3, height = 2)
