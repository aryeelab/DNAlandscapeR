library(diffloop)

load("/Users/lareauc/Desktop/Research/AryeeResearch/sarah-qual/output/valid_full-v1.0.rda")
ctcf <- rmchr(padGRanges(bedToGRanges("/Users/lareauc/Desktop/Research/AryeeResearch/sarah-qual/input/CTCF-np_peaks.narrowPeak"), pad = 1000))
h3k27ac <- rmchr(padGRanges(bedToGRanges("/Users/lareauc/Desktop/Research/AryeeResearch/sarah-qual/input/H3K27ac-np_peaks.narrowPeak"), pad = 1000))
promoter <- padGRanges(getHumanTSS(), pad = 1000)
vf_anno <- annotateLoops(valid_full, ctcf, h3k27ac, promoter)

# Very important to match data object name with the prefix of the .rda

Naive_1_SMC1 <- vf_anno[,1]
Naive_2_SMC1 <- vf_anno[,2]
Primed_1_SMC1 <- vf_anno[,3]
Primed_2_SMC1 <- vf_anno[,4]
Jurkat_1_SMC1 <- vf_anno[,5]
Jurkat_2_SMC1 <- vf_anno[,6]

save(Naive_1_SMC1, file="../data/loops/Naive_1_SMC1.rda")
save(Naive_2_SMC1, file="../data/loops/Naive_2_SMC1.rda")
save(Primed_1_SMC1, file="../data/loops/Primed_1_SMC1.rda")
save(Primed_2_SMC1, file="../data/loops/Primed_2_SMC1.rda")
save(Jurkat_1_SMC1, file="../data/loops/Jurkat_1_SMC1.rda")
save(Jurkat_2_SMC1, file="../data/loops/Jurkat_2_SMC1.rda")