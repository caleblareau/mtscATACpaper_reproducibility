library(data.table)
library(BuenColors)

d <- fread("../output/barcode_qc/CD34_H8.barcode_qc.tsv")

qplot(log10(d$depth),d$FRIP)

qplot(log10(d$depth), log10(d$mtDNAcoverage))
