library(rtracklayer)
library(GenomeInfoDb)
si <- Seqinfo(genome="mm9")

#wigToBigWig("/Users/lareauc/Downloads/GSM1397343_06182009_42A6WAAXX_B6.wig.gz", seqinfo = si, dest = "../data/mouse/tracks/ESC-H3K37me3.bw", clip = TRUE)