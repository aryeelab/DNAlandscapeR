library(diffloop)

## Read one sample at a time
promoter <- padGRanges(getHumanTSS(), pad = 1000)

x <- loopsMake.mango(beddir = "/Users/lareauc/Desktop/Research/AryeeResearch/processed_chiapet/human/POL2",
                       samples = "gm12878-pol2_tang", ext = "all")
enhancer <- padGRanges(bedToGRanges("/Users/lareauc/Desktop/GSM733771_hg19_wgEncodeBroadHistoneGm12878H3k27acStdPk.broadPeak"),  pad = 1000)
y <- annotateLoops(x, enhancer = rmchr(enhancer), promoter = promoter)
saveRDS(y, "/Users/lareauc/Desktop/GM12878-ChIA-Pet-POL2.rds")
