# Author: Caleb Lareau
# Group:  Aryee Lab

# Use: Creates list of chromosomal matrices
# of Hi-C data in .rds object for DNAlandscapeR

library(readr)
library(GenomicRanges)
library(reshape2)
library(tools)

matrix.file <- "gm12878_1000000.matrix"
bed.file <- "gm12878_1000000_abs.bed"
out.pre <-  basename(file_path_sans_ext(matrix.file))
compress = FALSE

bed.GRanges <- GRanges(data.frame(read_tsv(bed.file, col_names = c("chr", "start", "stop", "region"))))
dat.long <- read_tsv(matrix.file, col_names = c("idx1", "idx2", "region"))

core_chrom <- paste("chr", seq(1,22,1), sep = "")
core_chrom <- c(core_chrom, "chrX")

dir.create(out.pre)
lapply(core_chrom, function(chr){
    cur.chrom <- bed.GRanges[seqnames(bed.GRanges) == chr]
    vals <- mcols(cur.chrom)$region
    dat.chrom <- dat.long[dat.long$idx1 %in% vals & dat.long$idx2 %in% vals,]
    mat.chrom <- dcast(data = dat.chrom, formula = idx2 ~ idx1, value.var = "region", fill = NA)
    row.names <- as.character(mat.chrom[,1])
    mat.chrom <- mat.chrom[,-1]

    starts <- start(cur.chrom)
    names(starts) <- as.character(vals)

    rownames(mat.chrom) <- unname(starts[row.names])
    colnames(mat.chrom) <- unname(starts[colnames(mat.chrom)])
    saveRDS(as.matrix(mat.chrom), file = paste(out.pre, "/", out.pre, "-", chr, ".rds", sep=""))
    print(chr)
})

if(compress) {tar(paste(out.pre,".tgz",sep=""),out.pre,compression='gzip'); unlink(out.pre, recursive = TRUE)}