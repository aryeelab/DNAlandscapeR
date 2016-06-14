library(biomaRt)

## HUMAN ##

core_chrom <- seq(1,22,1)
core_chrom <- c(core_chrom, "X")

chr <- core_chrom

# Human protein coding genes
vals = list(chr, "protein_coding")
mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "grch37.ensembl.org", 
    path = "/biomart/martservice", dataset = "hsapiens_gene_ensembl")
geneinfo = getBM(attributes = c("chromosome_name", "start_position", 
    "end_position", "external_gene_name", "strand"), filters = c("chromosome_name", 
    "biotype"), values = vals, mart = mart)
colnames(geneinfo) <- c("chrom", "start", "stop", "gene", "strand")
geneinfo$score <- "."
geneinfo <- geneinfo[,c(1,2,3,4,6,5)]
save(geneinfo, file="data/GenomeAnnotation/hg19/geneinfo.rda")

# Human all exons
vals = list(chr)
mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org",
    path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
geneinfo = getBM(attributes = c("chromosome_name", "exon_chrom_start", "exon_chrom_end", "external_gene_name", 
    "strand"), mart = mart, values = vals)
colnames(geneinfo) <- c("chrom", "start", "stop", "gene", "strand")
geneinfo <- geneinfo[geneinfo$chrom %in% chr,]
geneinfo$score <- "."
geneinfo <- geneinfo[,c(1,2,3,4,6,5)]
row.names(geneinfo) <- NULL
save(geneinfo, file="data/GenomeAnnotation/hg19/geneinfo-exon.rda")

## MOUSE ##

core_chrom <- seq(1,19,1)
core_chrom <- c(core_chrom, "X")

chr <- core_chrom

# Mouse protein coding
vals = list(chr, "protein_coding")
mart = useMart(host='may2012.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL',
             dataset = "mmusculus_gene_ensembl")
geneinfo = getBM(attributes = c("chromosome_name", "start_position", 
    "end_position", "external_gene_id", "strand"), filters = c("chromosome_name", 
    "biotype"), values = vals, mart = mart)
colnames(geneinfo) <- c("chrom", "start", "stop", "gene", "strand")
geneinfo$score <- "."
geneinfo <- geneinfo[,c(1,2,3,4,6,5)]
save(geneinfo, file="data/GenomeAnnotation/mm9/geneinfo.rda")

# Mouse all exons
vals = list(chr, "protein_coding")
mart = useMart(host='may2012.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL',
             dataset = "mmusculus_gene_ensembl")
geneinfo = getBM(attributes = c("chromosome_name", "exon_chrom_start", 
    "exon_chrom_end", "external_gene_id", "strand"), filters = c("chromosome_name"), values = vals, mart = mart)
colnames(geneinfo) <- c("chrom", "start", "stop", "gene", "strand")
geneinfo <- geneinfo[geneinfo$chrom %in% chr,]
geneinfo$score <- "."
geneinfo <- geneinfo[,c(1,2,3,4,6,5)]
row.names(geneinfo) <- NULL
save(geneinfo, file="data/GenomeAnnotation/mm9/geneinfo-exon.rda")

