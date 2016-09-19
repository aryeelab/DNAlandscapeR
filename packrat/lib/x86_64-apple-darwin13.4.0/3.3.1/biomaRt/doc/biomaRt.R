### R code from vignette source 'biomaRt.Rnw'

###################################################
### code chunk number 1: annotate
###################################################
## library("annotate")
options(width=120)


###################################################
### code chunk number 2: biomaRt
###################################################
library("biomaRt")
listMarts()


###################################################
### code chunk number 3: ensembl1
###################################################
ensembl=useMart("ensembl")


###################################################
### code chunk number 4: listDatasets
###################################################
listDatasets(ensembl)


###################################################
### code chunk number 5: ensembl2
###################################################
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")


###################################################
### code chunk number 6: filters
###################################################
filters = listFilters(ensembl)
filters[1:5,]


###################################################
### code chunk number 7: attributes
###################################################
attributes = listAttributes(ensembl)
attributes[1:5,]


###################################################
### code chunk number 8: biomaRt.Rnw:120-122 (eval = FALSE)
###################################################
## affyids=c("202763_at","209310_s_at","207500_at")
## getBM(attributes=c('affy_hg_u133_plus_2', 'entrezgene'), filters = 'affy_hg_u133_plus_2', values = affyids, mart = ensembl)


###################################################
### code chunk number 9: biomaRt.Rnw:140-143 (eval = FALSE)
###################################################
## affyids=c("202763_at","209310_s_at","207500_at")
## getBM(attributes=c('affy_hg_u133_plus_2', 'hgnc_symbol', 'chromosome_name','start_position','end_position', 'band'),
##  filters = 'affy_hg_u133_plus_2', values = affyids, mart = ensembl)


###################################################
### code chunk number 10: biomaRt.Rnw:158-161 (eval = FALSE)
###################################################
## entrez=c("673","837")
## goids = getBM(attributes=c('entrezgene','go_id'), filters='entrezgene', values=entrez, mart=ensembl)
## head(goids)


###################################################
### code chunk number 11: biomaRt.Rnw:197-199 (eval = FALSE)
###################################################
## refseqids = c("NM_005359","NM_000546")
## ipro = getBM(attributes=c("refseq_mrna","interpro","interpro_description"), filters="refseq_mrna",values=refseqids, mart=ensembl)


###################################################
### code chunk number 12: biomaRt.Rnw:221-223
###################################################
getBM(c('affy_hg_u133_plus_2','ensembl_gene_id'), filters = c('chromosome_name','start','end'),
 values=list(16,1100000,1250000), mart=ensembl)


###################################################
### code chunk number 13: biomaRt.Rnw:230-231 (eval = FALSE)
###################################################
## getBM(c('entrezgene','hgnc_symbol'), filters='go', values='GO:0004707', mart=ensembl)


###################################################
### code chunk number 14: biomaRt.Rnw:252-254 (eval = FALSE)
###################################################
## entrez=c("673","7157","837")
## getSequence(id = entrez, type="entrezgene",seqType="coding_gene_flank",upstream=100, mart=ensembl) 


###################################################
### code chunk number 15: biomaRt.Rnw:262-265 (eval = FALSE)
###################################################
## utr5 = getSequence(chromosome=3, start=185514033, end=185535839,
##                       type="entrezgene",seqType="5utr", mart=ensembl)
## utr5


###################################################
### code chunk number 16: biomaRt.Rnw:282-285 (eval = FALSE)
###################################################
## protein = getSequence(id=c(100, 5728),type="entrezgene",
##                         seqType="peptide", mart=ensembl)
## protein


###################################################
### code chunk number 17: biomaRt.Rnw:301-302 (eval = FALSE)
###################################################
## snpmart = useMart("snp", dataset="hsapiens_snp")


###################################################
### code chunk number 18: biomaRt.Rnw:309-310 (eval = FALSE)
###################################################
## getBM(c('refsnp_id','allele','chrom_start','chrom_strand'), filters = c('chr_name','chrom_start','chrom_end'), values = list(8,148350,148612), mart = snpmart)


###################################################
### code chunk number 19: biomaRt.Rnw:363-364
###################################################
listMarts(archive=TRUE)


###################################################
### code chunk number 20: biomaRt.Rnw:369-370 (eval = FALSE)
###################################################
## ensembl = useMart("ensembl_mart_46", dataset="hsapiens_gene_ensembl", archive = TRUE)


###################################################
### code chunk number 21: biomaRt.Rnw:381-384 (eval = FALSE)
###################################################
## listMarts(host='may2009.archive.ensembl.org')
## ensembl54=useMart(host='may2009.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL')
## ensembl54=useMart(host='may2009.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl')


###################################################
### code chunk number 22: biomaRt.Rnw:393-400 (eval = FALSE)
###################################################
## wormbase=useMart("WS220",dataset="wormbase_gene")
## listFilters(wormbase)
## listAttributes(wormbase)
## getBM(attributes = c("public_name","rnai","rnai_phenotype_phenotype_label"),
##                      filters="gene_name", values=c("unc-26","his-33"),
##                      mart=wormbase)
##      


###################################################
### code chunk number 23: biomaRt.Rnw:432-433
###################################################
filterType("with_affy_hg_u133_plus_2",ensembl)


###################################################
### code chunk number 24: biomaRt.Rnw:442-443
###################################################
filterOptions("biotype",ensembl)


###################################################
### code chunk number 25: biomaRt.Rnw:458-460
###################################################
pages = attributePages(ensembl)
pages


###################################################
### code chunk number 26: biomaRt.Rnw:467-468
###################################################
listAttributes(ensembl, page="feature_page")


###################################################
### code chunk number 27: columnsAndKeyTypes
###################################################
mart<-useMart(dataset="hsapiens_gene_ensembl",biomart='ensembl')
head(keytypes(mart), n=3)
head(columns(mart), n=3)


###################################################
### code chunk number 28: keys1
###################################################
k = keys(mart, keytype="chromosome_name")
head(k, n=3)


###################################################
### code chunk number 29: keys2
###################################################
k = keys(mart, keytype="chromosome_name", pattern="LRG")
head(k, n=3)


###################################################
### code chunk number 30: select
###################################################
affy=c("202763_at","209310_s_at","207500_at")
select(mart, keys=affy, columns=c('affy_hg_u133_plus_2','entrezgene'),
  keytype='affy_hg_u133_plus_2')


###################################################
### code chunk number 31: biomaRt.Rnw:554-556
###################################################
sessionInfo()
warnings()


