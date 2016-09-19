## ----pre,echo=FALSE,results='hide'---------------------------------------
library(knitr)
opts_chunk$set(warning=FALSE,message=FALSE,cache=TRUE)

## ----loadLibrary---------------------------------------------------------
library(GEOquery)

## ------------------------------------------------------------------------
# If you have network access, the more typical way to do this
# would be to use this:
# gds <- getGEO("GDS507")
gds <- getGEO(filename=system.file("extdata/GDS507.soft.gz",package="GEOquery"))

## ------------------------------------------------------------------------
# If you have network access, the more typical way to do this
# would be to use this:
# gds <- getGEO("GSM11805")
gsm <- getGEO(filename=system.file("extdata/GSM11805.txt.gz",package="GEOquery"))

## ------------------------------------------------------------------------
# Look at gsm metadata:
head(Meta(gsm))
# Look at data associated with the GSM:
# but restrict to only first 5 rows, for brevity
Table(gsm)[1:5,]
# Look at Column descriptions:
Columns(gsm)

## ------------------------------------------------------------------------
Columns(gds)[,1:3]

## ------------------------------------------------------------------------
# Again, with good network access, one would do:
# gse <- getGEO("GSE781",GSEMatrix=FALSE)
gse <- getGEO(filename=system.file("extdata/GSE781_family.soft.gz",package="GEOquery"))
head(Meta(gse))
# names of all the GSM objects contained in the GSE
names(GSMList(gse))
# and get the first GSM object on the list
GSMList(gse)[[1]]
# and the names of the GPLs represented
names(GPLList(gse))

## ------------------------------------------------------------------------
# Note that GSEMatrix=TRUE is the default
gse2553 <- getGEO('GSE2553',GSEMatrix=TRUE)
show(gse2553)
show(pData(phenoData(gse2553[[1]]))[1:5,c(1,6,8)])

## ------------------------------------------------------------------------
eset <- GDS2eSet(gds,do.log2=TRUE)

## ------------------------------------------------------------------------
eset
pData(eset)[,1:3]

## ------------------------------------------------------------------------
#get the platform from the GDS metadata
Meta(gds)$platform
#So use this information in a call to getGEO
gpl <- getGEO(filename=system.file("extdata/GPL97.annot.gz",package="GEOquery"))

## ------------------------------------------------------------------------
MA <- GDS2MA(gds,GPL=gpl)
class(MA)

## ------------------------------------------------------------------------
gsmplatforms <- lapply(GSMList(gse),function(x) {Meta(x)$platform})
head(gsmplatforms)

## ------------------------------------------------------------------------
gsmlist = Filter(function(gsm) {Meta(gsm)$platform=='GPL96'},GSMList(gse))
length(gsmlist)

## ------------------------------------------------------------------------
Table(gsmlist[[1]])[1:5,]
# and get the column descriptions
Columns(gsmlist[[1]])[1:5,]

## ------------------------------------------------------------------------
# get the probeset ordering
probesets <- Table(GPLList(gse)[[1]])$ID
# make the data matrix from the VALUE columns from each GSM
# being careful to match the order of the probesets in the platform
# with those in the GSMs
data.matrix <- do.call('cbind',lapply(gsmlist,function(x) 
                                      {tab <- Table(x)
                                       mymatch <- match(probesets,tab$ID_REF)
                                       return(tab$VALUE[mymatch])
                                     }))
data.matrix <- apply(data.matrix,2,function(x) {as.numeric(as.character(x))})
data.matrix <- log2(data.matrix)
data.matrix[1:5,]

## ------------------------------------------------------------------------
require(Biobase)
# go through the necessary steps to make a compliant ExpressionSet
rownames(data.matrix) <- probesets
colnames(data.matrix) <- names(gsmlist)
pdata <- data.frame(samples=names(gsmlist))
rownames(pdata) <- names(gsmlist)
pheno <- as(pdata,"AnnotatedDataFrame")
eset2 <- new('ExpressionSet',exprs=data.matrix,phenoData=pheno)
eset2

## ------------------------------------------------------------------------
gpl97 <- getGEO('GPL97')
Meta(gpl97)$title
head(Meta(gpl97)$series_id)
length(Meta(gpl97)$series_id)
head(Meta(gpl97)$sample_id)
length(Meta(gpl97)$sample_id)

## ------------------------------------------------------------------------
gsmids <- Meta(gpl97)$sample_id
gsmlist <- sapply(gsmids[1:5],getGEO)
names(gsmlist)

## ----citation------------------------------------------------------------
citation("GEOquery")

## ----eval=FALSE----------------------------------------------------------
#  bug.report(package='GEOquery')

## ----echo=FALSE----------------------------------------------------------
sessionInfo()

