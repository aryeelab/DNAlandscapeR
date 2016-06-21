# Collate Hi-C data and run association testing on it

library(GenomicRanges)
library(foreach)
library(Matrix)
library(utils)
library(matrixStats)
library(data.table)
library(edgeR)
library(reshape2)
options(scipen=999)


samples <- c("BT142-rep1", "BT142-rep2", "GBM4-rep1", "GBM4-rep2")
groups <- c("BT142", "BT142", "GBM4", "GBM4")
res <- "500000" 
pre <- paste0("GBM4-BT142diff_", res)
compress = TRUE


buildInteractionsCountsMatrix <- function(samples, res, chr){
    dat <- foreach(sample = samples) %do% readRDS( paste0(sample, "/", sample, "_", res, "/", sample, "_", res, "-", chr, ".rds"))
    nonZero <- summary(Reduce("*", dat))
    counts <- sapply(dat, function(m){ as.matrix(m[cbind(nonZero$i, nonZero$j)])})
    names <- colnames(dat[[1]])
    joined <- data.frame(names[nonZero$i], names[nonZero$j], counts, chr, stringsAsFactors = FALSE)
    colnames(joined) <- c("region1", "region2", samples, "chrom")
    
    return(joined)
}

core_chrom <- paste("chr", seq(1, 22, 1), sep = "")
core_chrom <- c(core_chrom, "chrX")
dist <- c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,
          135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,
          59128983,63025520,48129895,51304566,155270560) # hg19 chromosome distances
names(dist) <- core_chrom

dat <- lapply(core_chrom, function(chr) { buildInteractionsCountsMatrix(samples, res, chr)})
names(dat) <- core_chrom
counts <- data.matrix(Reduce(rbind, dat)[,c(3:(2+ length(samples)))])

# Run edgeR
z <- DGEList(counts = counts, group = groups)
design <- model.matrix(~groups)
z <- calcNormFactors(z)
yy <- estimateDisp(z, design)
fit <- glmQLFit(yy, design, robust = TRUE)
qlf <- glmQLFTest(fit, coef = 2)
results <- as.data.frame(topTags(qlf, n = nrow(counts), sort.by = "none"))

# Print table of top results
int_res <- cbind(Reduce(rbind, dat)[,c(length(samples) + 3, 2,1)], results, counts)
int_top_res <- int_res[with(int_res, order(FDR)), ]; int_top_res <- int_top_res[int_top_res$FDR < 0.2, ]
#int_top_res <- head(int_res[with(int_res, order(FDR)), ], 100)
write.table(int_top_res, sep = "\t", row.names = FALSE, quote = FALSE, file = paste0(pre, ".tsv"))

# Re-create sparse matrices per chromosome
dir.create(pre)
n_row <- 0
lapply(core_chrom, function(chr) {
    
    chrom.df <- dat[[chr]]
    dat.chrom <- data.matrix(cbind(as.numeric(chrom.df[,2]), as.numeric(chrom.df[,1]), results[c((n_row + 1):(dim(chrom.df)[1])), ]$logFC))
    colnames(dat.chrom) <- c("idx1", "idx2", "region")
    n_row <- n_row + dim(chrom.df)[1]
    
    # Create zeroes matrix to handle missing data
    bins <- seq(0, as.numeric(dist[chr]), as.numeric(res))
    zeros.long <- cbind(t(combn(bins, 2)), 0)
    
    zeros.long <- rbind(zeros.long, cbind(bins, bins, 0))
    colnames(zeros.long) <- c("idx1", "idx2", "region")
    
    mat.chrom <- dcast(data = data.frame(rbind(dat.chrom, zeros.long)), formula = idx2 ~ idx1,
                       value.var = "region", fill = 0, fun.aggregate = sum)
    row.names(mat.chrom) <- as.character(mat.chrom[, 1])
    mat.chrom <- mat.chrom[, -1]
    
    i <- dim(mat.chrom)[1]
    j <- dim(mat.chrom)[2]
    for (k in 1:i) mat.chrom[(k + 41):j, k] <- 0
    mat.chrom <- mat.chrom[1:i, 1:j]
    
    print(dim(mat.chrom))
    saveRDS(Matrix(as.matrix(mat.chrom)), file = paste(pre, "/", pre, "-", chr, ".rds", sep = ""))
})

if (compress) { tar(paste(pre, ".tgz", sep = ""), pre, compression = 'gzip'); unlink(pre, recursive = TRUE) }
