library(data.table)
library(dplyr)

compute_nuclear_AFs <- function(hsdf){
  lef_ref <- hsdf %>% filter(V2 == 109084804 & V3 == "A") %>% dim
  lef_alt <- hsdf %>% filter(V2 == 109084804 & V3 == "C") %>% dim
  hcst_ref <- hsdf %>% filter(V2 == 36394730 & V3 == "G") %>% dim
  hcst_alt <- hsdf %>% filter(V2 == 36394730 & V3 == "A") %>% dim
  
  data.frame(
    LEF1 = lef_alt[1]/(lef_alt[1] + lef_ref[1])*100,
    HCST = hcst_alt[1]/(hcst_alt[1] + hcst_ref[1])*100
  )
  
}

cd19pos <- data.frame(fread("../scRNAseq_data/nuclear_mutations_pt2/hisnipper_CD19pos.tsv"))
cd19neg <- data.frame(fread("../scRNAseq_data/nuclear_mutations_pt2/hisnipper_CD19neg.tsv"))

# compute allele frequencies for supplemental figure
compute_nuclear_AFs(cd19pos)
compute_nuclear_AFs(cd19neg)
