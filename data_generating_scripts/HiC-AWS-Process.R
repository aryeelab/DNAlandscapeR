# Author: Caleb Lareau
# Group:  Aryee Lab

# Use: Creates list of chromosomal matrices
# of Hi-C data in .rds object for DNAlandscapeR

# v2.0-- Creates full (coordinates) sparse matrix

library(readr)
library(GenomicRanges)
library(reshape2)
library(Matrix)
library(foreach)

b <- "/PHShome/ma695/work/projects/gbm_topology/output/hicpro/johnstone_gbm/hic_results/matrix/"
s_all <- c("BT142-rep1", "BT142-rep2", "GBM4-rep1", "GBM4-rep2")
c_all <- c("20000", "40000", "150000", "500000", "1000000")

foreach( s = s_all ) %do% { foreach( c = c_all) %do% {
matrix.file <- paste0(b,s, "/iced/", c, "/",s,"_", c, "_iced.matrix")
bed.file <- paste0(b,s, "/raw/", c, "/",s,"_", c, "_abs.bed")
out.pre <- paste0(s, "_", c)
print(out.pre)
compress = TRUE

bed.GRanges <- GRanges(data.frame(read_tsv(bed.file, col_names = c("chr", "start", "stop", "region"))))
dat.long <- read_tsv(matrix.file, col_names = c("idx1", "idx2", "region"))

core_chrom <- paste("chr", seq(1, 22, 1), sep = "")
core_chrom <- c(core_chrom, "chrX")
dist <- c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,
          135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,
          59128983,63025520,48129895,51304566,155270560) # hg19 chromosome distances
names(dist) <- core_chrom

dir.create(out.pre)
lapply(core_chrom, function(chr) {
  cur.chrom <- bed.GRanges[seqnames(bed.GRanges) == chr]
  vals <-  mcols(cur.chrom)$region 
  dat.chrom <- dat.long[dat.long$idx1 %in% vals & dat.long$idx2 %in% vals, ]
  starts <- start(cur.chrom)
  names(starts) <- vals
  dat.chrom$idx1 <- starts[as.character(dat.chrom$idx1)]
  dat.chrom$idx2 <- starts[as.character(dat.chrom$idx2)]
  
  # Create zeroes matrix to handle missing data
  bins <- seq(0, as.numeric(dist[chr]), as.numeric(c))
  zeros.long <- cbind(t(combn(bins, 2)), 0)
  zeros.long <- rbind(zeros.long, cbind(bins, bins, 0))
  colnames(zeros.long) <- c("idx1", "idx2", "region")
  
  mat.chrom <-dcast(data = rbind(dat.chrom, zeros.long), formula = idx2 ~ idx1,
                    value.var = "region", fill = 0, fun.aggregate = sum)
  row.names(mat.chrom) <- as.character(mat.chrom[, 1])
  mat.chrom <- mat.chrom[, -1]
  
  i <- dim(mat.chrom)[1]
  j <- dim(mat.chrom)[2]
  for (k in 1:i) mat.chrom[(k + 41):j, k] <- 0
  mat.chrom <- mat.chrom[1:i, 1:j]
  
  saveRDS(Matrix(as.matrix(mat.chrom)), file = paste(out.pre, "/", out.pre, "-", chr, ".rds", sep = ""))
  print(chr)
})

if(compress) {tar(paste0(out.pre,".tgz"),out.pre,compression='gzip'); unlink(out.pre, recursive = TRUE)}
}}