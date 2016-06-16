# Author: Caleb Lareau
# Group:  Aryee Lab

# Use: Creates list of sparse matrices of Hi-C data indexed by chromosome in a
# .rds object for visualization in DNAlandscapeR

library(readr)
library(GenomicRanges)
library(reshape2)
library(tools)
library(Matrix)

matrix.file <- "gm12878_1000000.matrix"
bed.file <- "gm12878_1000000_abs.bed"
out.prefix <-  "gm12878_1000000"

bed.GRanges <- GRanges(data.frame(read_tsv(bed.file, col_names = c("chr", "start", "stop", "region"))))
dat.long <- read_tsv(matrix.file, col_names = c("idx1", "idx2", "region"))

core_chrom <- paste0("chr", seq(1,22,1))
core_chrom <- c(core_chrom, "chrX")

dat <- lapply(core_chrom, function(chr){
  cur.chrom <- bed.GRanges[seqnames(bed.GRanges) == chr]
  vals <- mcols(cur.chrom)$region
  dat.chrom <- dat.long[dat.long$idx1 %in% vals & dat.long$idx2 %in% vals,]
  mat.chrom <- dcast(data = dat.chrom, formula = idx2 ~ idx1, value.var = "region", fill = 0)
  row.names <- as.character(mat.chrom[,1])
  mat.chrom <- mat.chrom[,-1]
  
  starts <- start(cur.chrom)
  names(starts) <- as.character(vals)
  rownames(mat.chrom) <- unname(starts[row.names])
  colnames(mat.chrom) <- unname(starts[colnames(mat.chrom)])
  
  i <- dim(mat.chrom)[1]; j <- dim(mat.chrom)[2]
  for(k in 1:i) mat.chrom[(k+41):j, k] <- 0
  mat.chrom <- mat.chrom[1:i, 1:j]
  Matrix(as.matrix(mat.chrom))
})

names(dat) <- core_chrom
saveRDS(dat, file = paste0(out.prefix, ".rds"))
