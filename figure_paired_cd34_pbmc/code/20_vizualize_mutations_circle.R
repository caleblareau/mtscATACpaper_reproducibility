library(SummarizedExperiment)
library(BuenColors)
library(reshape2)
library(stringr)

"%ni%" <- Negate("%in%")

# Import per-sample calls
pbmc <- readRDS("../output/filteredpbmcs_mgatk_calls.rds")
cd34 <- readRDS("../output/filteredCD34_mgatk_calls.rds")

plot_df <- data.frame(
  pos = rowData(cd34)$position,
  nucleotide = rowData(cd34)$nucleotide,
  newAF = (7474*rowData(cd34)$mean + 8591*rowData(pbmc)$mean)/(8591 + 7474)
)
plot_df$color_code = plot_df$nucleotide %in% c("G>A", "A>G", "C>T", "T>C")

# Modify Stephen Turner's code
addgenelabel <- function(bp,gene) { gene <- ifelse(bp < 577,gene <- "Control-Region", ifelse(bp < 648,gene <- "tRNA", ifelse(bp < 1602,gene <- "rRNA", ifelse(bp < 1671,gene <- "tRNA", ifelse(bp < 3230,gene <- "rRNA", ifelse(bp < 3305,gene <- "tRNA", ifelse(bp < 3307,gene <- "Non-Coding", ifelse(bp < 4263,gene<- "ND1", ifelse(bp < 4332,gene <- "tRNA", ifelse(bp < 4401,gene <- "tRNA", ifelse(bp < 4402,gene <- "Non-Coding", ifelse(bp < 4470,gene <- "tRNA", ifelse(bp < 5512,gene <- "ND2", ifelse(bp < 5580,gene <- "tRNA", ifelse(bp < 5587,gene <- "Non-Coding", ifelse(bp < 5656,gene <- "tRNA", ifelse(bp < 5657,gene <- "Non-Coding", ifelse(bp < 5730,gene <- "tRNA", ifelse(bp < 5826,gene <- "tRNA", ifelse(bp < 5892,gene <- "tRNA", ifelse(bp < 5904,gene <- "Non-Coding", ifelse(bp < 7446,gene <- "CO1", ifelse(bp < 7515,gene <- "tRNA", ifelse(bp < 7518,gene <- "Non-Coding", ifelse(bp < 7586,gene <- "tRNA", ifelse(bp < 8270,gene <- "CO2", ifelse(bp < 8295,gene <- "Non-Coding", ifelse(bp < 8365,gene <- "tRNA", ifelse(bp < 8366,gene <- "Non-Coding", ifelse(bp < 8573,gene <- "ATP8", ifelse(bp < 9208,gene <- "ATP6", ifelse(bp < 9991,gene <- "CO3", ifelse(bp < 10059,gene <- "tRNA", ifelse(bp < 10405,gene <- "ND3", ifelse(bp < 10470,gene <- "tRNA", ifelse(bp < 10767,gene <- "ND4L", ifelse(bp < 12138,gene <- "ND4", ifelse(bp < 12207,gene <- "tRNA", ifelse(bp < 12266,gene <- "tRNA", ifelse(bp < 12337,gene <- "tRNA", ifelse(bp < 14149,gene <- "ND5", ifelse(bp < 14674,gene <- "ND6", ifelse(bp < 14743,gene <- "tRNA", ifelse(bp < 14747,gene <- "Non-Coding", ifelse(bp < 15888,gene <- "CYB", ifelse(bp < 15954,gene <- "tRNA", ifelse(bp < 15956,gene <- "Non-Coding", ifelse(bp < 16024,gene <- "tRNA", ifelse(bp < 17000,gene <- "Control-Region") ))))))))))))))))))) ))))))))))))))))))) ))))))))) ) }
visibleboundaries <- c(576,1601,3229,4262,5511,7445,8269,9207,9990,10404,10766,12137,14148,14673,15887)

lines <- data.frame(x = seq(0,16569,by=1),y = 0)
lines$gene <- addgenelabel(lines$x,lines$gene)

colours <- c("Control-Region" = "lightblue4", "tRNA" = "black", "rRNA" = "mediumaquamarine", "Non-Coding" ="sienna4", "ND1" = "magenta", "ND2" = "mediumblue", "CO1" = "olivedrab", "CO2" = "orange2", "ATP8" = "orchid4", "ATP6" = "red3",
             "CO3" = "royalblue2", "ND3" = "palegreen4", "ND4L" = "grey0", "ND4" = "pink4", "ND5" = "yellow4", "ND6" = "steelblue4", "CYB" ="tan", "black")

gene_points <- data.frame(pos = 1:16569, gene = "" )
gene_points$color_code  <- addgenelabel(gene_points$pos, gene_points$gene)
gene_points$newAF <- 0.0001

xxtheme <-   theme(
  axis.line = element_blank(),
  #axis.ticks.y = element_blank(),        ## <- this line
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  panel.border = element_blank()) 

thin <- c(rep(FALSE,10),TRUE,rep(FALSE,10))

special <- c("TRUE" = "black", "FALSE" = "dodgerblue2")
p1 <- ggplot(plot_df, aes(x = pos, y = newAF, color = color_code)) + geom_point(size = 0.3, color = "black") + 
  pretty_plot(fontsize = 6)  + geom_point(data = gene_points[thin,], size = 0.2)+
  scale_color_manual(values = c(special, colours)) +
  coord_polar(direction = 1) + labs(x = "", y = "Heteroplasmy") + 
  scale_y_log10(breaks=c(0.00001,0.0001,0.001, 0.01), limits = c(0.00001, 0.01)) + theme(legend.position = "none") + xxtheme 

cowplot::ggsave2(p1, filename = "../plots/circos_all_mutations.pdf", height = 3, width = 3)

