#' Process Jurkat 450k methylation

library(minfi)

dir <- "/Users/lareauc/Downloads/jurkat-450k"
tb <- file.path(dir, c("GSM999367_hg19_wgEncodeHaibMethyl450JurkatSitesRep1"))
rgset <- read.metharray(tb, verbose = TRUE)
mset <- preprocessIllumina(rgset)
mset <- mapToGenome(mset)
loc <- data.frame(granges(mset))
loc$val <- getBeta(mset)
loc <- loc[c(-4,-5)]
colnames <- c("chrom", "start", "end", "val") 
write.table(loc, file = "jurkat450k.bedgraph")