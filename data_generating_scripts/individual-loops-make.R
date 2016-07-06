library(diffloop)

## Read one sample at a time
promoter <- padGRanges(getHumanTSS(), pad = 1000)


x <- loopsMake.mango(beddir = "/Users/lareauc/Desktop/Research/AryeeResearch/processed_chiapet/human/POL2",
                       samples = "mcf7-pol2_c", ext = "all")
enhancer <- padGRanges(bedToGRanges("/Users/lareauc/Downloads/GSM945854_hg19_wgEncodeSydhHistoneMcf7H3k27acUcdPk.narrowPeak"),  pad = 1000)
y <- annotateLoops(x, enhancer = rmchr(enhancer), promoter = promoter)
saveRDS(y, "/Users/lareauc/Desktop/MCF7-ChIA-Pet-POL2.rds")
