library(data.table)
library(diffloop)
"%ni%" <- Negate("%in%")

# Import for {f}ixed or {r}egular for {t}f1 or {m}errf
f_t <- bedToGRanges("m2_fixed_TF1_peaks.narrowPeak.gz")
f_m <- bedToGRanges("m2_fixed_MERRF_peaks.narrowPeak.gz")

r_t <- bedToGRanges("m2_regular_TF1_peaks.narrowPeak.gz")
r_m <- bedToGRanges("m2_regular_MERRF_peaks.narrowPeak.gz")

# Overlap regulars to find CTS speaks
ov_reg <- findOverlaps(r_t, r_m)
cts <- c(r_t[1:length(r_t) %ni% queryHits(ov_reg)], 
         r_m[1:length(r_m) %ni% subjectHits(ov_reg)])

# Determine the number of CTS peaks found by the regular protocol
length(unique(subjectHits(findOverlaps(c(f_t, f_m), cts))))
length(cts)

# Get the proportion to report
length(unique(subjectHits(findOverlaps(c(f_t, f_m), cts))))/length(cts) *100
