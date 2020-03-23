library(data.table)
library(rtracklayer)
"%ni%" <- Negate("%in%")

cdf <- readRDS("bcell_substructure/output/PT1_clone_definition.rds")
frags <- fread("CLL_CD19pos_CR-mtMask/outs/fragments.tsv.gz", header = FALSE) %>% data.frame()

clusters <- as.character(unique(cdf$cluster_id))

sapply(clusters, function(cluster){
  
  possible_ids <- cdf %>% filter(cluster_id == cluster) %>% pull(cell_id) %>% as.character()
  cluster_gr <- frags %>% filter(V4 %in% possible_ids) %>%
    setnames(c("chr", "start", "end", "V4", "PCRn")) %>%
    makeGRangesFromDataFrame()
  
  reads_coverage <- coverage(cluster_gr)/length(cluster_gr)*1000000
  export.bw(reads_coverage, con = paste0("bigwigs/CLL-Patient1-Bcell-", as.character(cluster), ".bw"))
  cluster
}) -> bulk2

x2 = readRDS("bcell_substructure/output/PT1_X2.rds")
df = x2 %>% mutate(region = paste0(chr, ":", start-1000, "-", end +1000)) %>% arrange(desc(obs_x2))

clonal_bc <- cdf %>% filter(cluster_id %ni% c("0", "1", "4")) %>% pull(cell_id) %>% as.character()
cluster_gr1 <- frags %>% filter(V4 %in% clonal_bc) %>%
  setnames(c("chr", "start", "end", "V4", "PCRn")) %>%
  makeGRangesFromDataFrame()

reads_coverage <- coverage(cluster_gr1)/length(cluster_gr1)*1000000
export.bw(reads_coverage, con = paste0("bigwigs/CLL-Patient1-Bcell-", as.character("other-clonal"), ".bw"))