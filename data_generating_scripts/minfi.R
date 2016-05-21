#' Process 450k methylation

library(minfi)

dir <- "/Users/lareauc/Desktop/Research/AryeeResearch/raw.dat/mcf7-450k"
tb <- file.path(dir, c("GSM999373_hg19_wgEncodeHaibMethyl450Mcf7SitesRep1"))
rgset <- read.metharray(tb, verbose = TRUE)
mset <- preprocessIllumina(rgset)
mset <- mapToGenome(mset)
loc <- data.frame(granges(mset))
loc$val <- getBeta(mset)
loc <- loc[c(-4,-5)]
colnames <- c("chrom", "start", "end", "val") 
write.table(loc, file = "MCF7-450k.bedgraph", row.names = FALSE, quote = FALSE)