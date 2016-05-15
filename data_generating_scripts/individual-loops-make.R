library(diffloop)

load("/Users/lareauc/Desktop/Research/AryeeResearch/sarah-qual/output/valid_full-v1.0.rda")
ctcf <- rmchr(padGRanges(bedToGRanges("/Users/lareauc/Desktop/Research/AryeeResearch/sarah-qual/input/CTCF-np_peaks.narrowPeak"), pad = 1000))
h3k27ac <- rmchr(padGRanges(bedToGRanges("/Users/lareauc/Desktop/Research/AryeeResearch/sarah-qual/input/H3K27ac-np_peaks.narrowPeak"), pad = 1000))
promoter <- padGRanges(getHumanTSS(), pad = 1000)
vf_anno <- annotateLoops(valid_full, ctcf, h3k27ac, promoter)
vf_anno@colData$groups <- c("esc", "esc", "esc", "esc", "jurkat", "jurkat")
vf_res <- quickAssoc(vf_anno)
esc_anno <- vf_anno[,c(1,2,3,4)]
esc_anno@colData$groups <- c("naive", "naive", "primed", "primed")
esc_res <- quickAssoc(esc_anno)

Naive_1_SMC1 <- vf_anno[,1]
Naive_2_SMC1 <- vf_anno[,2]
Primed_1_SMC1 <- vf_anno[,3]
Primed_2_SMC1 <- vf_anno[,4]
Jurkat_1_SMC1 <- vf_anno[,5]
Jurkat_2_SMC1 <- vf_anno[,6]

saveRDS(Naive_1_SMC1, file="../data/loops/Naive_1_SMC1.rds")
saveRDS(Naive_2_SMC1, file="../data/loops/Naive_2_SMC1.rds")
saveRDS(Primed_1_SMC1, file="../data/loops/Primed_1_SMC1.rds")
saveRDS(Primed_2_SMC1, file="../data/loops/Primed_2_SMC1.rds")
saveRDS(Jurkat_1_SMC1, file="../data/loops/Jurkat_1_SMC1.rds")
saveRDS(Jurkat_2_SMC1, file="../data/loops/Jurkat_2_SMC1.rds")