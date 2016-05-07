library(rtracklayer)
library(GenomeInfoDb)
si <- Seqinfo(genome="hg19")

#wigToBigWig("/Users/lareauc/Downloads/GSM1705251_20150206_3642.4518.rpm.WIG", seqinfo = si, dest = "../data/Naive-R1-H3K27ac.bw", clip = TRUE)
#wigToBigWig("/Users/lareauc/Downloads/GSM1705260_20150206_3646.4478.rpm.WIG.gz", seqinfo = si, dest = "data/Primed-R1-H3K27ac.bw", clip = TRUE)
#wigToBigWig("/Users/lareauc/Downloads/GSM1697882_20141230_3391.1796.rpm.WIG.gz", seqinfo = si, dest = "data/Jurkat-H3K27ac.bw", clip = TRUE)
wigToBigWig("/Users/lareauc/Downloads/Naive_RNAseq.wig", seqinfo = si, dest = "data/Naive_RNAseq.bw")
a<-import("/Users/lareauc/Downloads/Naive_RNAseq.wig", format = "WIG", seqinfo = si)
gr1c <- as(as(slice(coverage(a), lower=1L, upper=1L), "GRanges"), "UCSCData")
export(head(a[unique(findOverlaps(a,gr1c)@from)]), format = "BigWig",  "data/Naive_RNAseq.bw")


##### METHYLATION DATA FROM SERVER ####
# load("bsseq-npdat.rda")
# library(GenomeInfoDb)
# library(rtracklayer)
# library(diffloop)
# si <- Seqinfo(genome="hg19")
# 
# dnc$score <- dnc$meth/(dnc$meth + dnc$unmeth)
# dnc$stop <- dnc$start
# dnc <- dnc[,c(-3,-4)]
# dnc <- dnc[dnc$chr!="MT",]
# u.dnc <- as(addchr(GRanges(dnc)), "UCSCData")
# u.dnc@seqinfo <- si
# export(u.dnc, format = "BigWig",  "Naive_methylation.bw")
# 
# dpc$score <- dpc$meth/(dpc$meth + dpc$unmeth)
# dpc$stop <- dpc$start
# dpc <- dpc[,c(-3,-4)]
# dpc <- dpc[dpc$chr!="MT",]
# u.dpc <- as(addchr(GRanges(dpc)), "UCSCData")
# u.dpc@seqinfo <- si
# export(u.dpc, format = "BigWig",  "Primed_methylation.bw")
