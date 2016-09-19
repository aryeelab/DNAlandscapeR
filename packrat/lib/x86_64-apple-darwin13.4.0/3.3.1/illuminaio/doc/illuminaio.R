### R code from vignette source 'illuminaio.Rnw'

###################################################
### code chunk number 1: style
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: loadData
###################################################
library(illuminaio)
library(IlluminaDataTestFiles)
idatFile <- system.file("extdata", "idat", "4343238080_A_Grn.idat",
                        package = "IlluminaDataTestFiles")
idat <- readIDAT(idatFile)                                 


###################################################
### code chunk number 3: exploreData
###################################################
names(idat)


###################################################
### code chunk number 4: printQuants
###################################################
idatData <- idat$Quants    
head(idatData)


###################################################
### code chunk number 5: readGenotyping
###################################################
genotypeIdatFile <- system.file("extdata", "idat", "5723646052_R02C02_Grn.idat",
                        package = "IlluminaDataTestFiles")
genotypeIdat <- readIDAT(genotypeIdatFile)
names(genotypeIdat)


###################################################
### code chunk number 6: printGenotypingQuants
###################################################
head(genotypeIdat$Quants)


###################################################
### code chunk number 7: ImportGenomeStudio
###################################################
gsFile <- system.file("extdata", "gs", "4343238080_A_ProbeSummary.txt.gz",
                      package = "IlluminaDataTestFiles")
gStudio <- read.delim(gsFile, sep = "\t", header = TRUE)
idatData <- idatData[which(idatData[,"CodesBinData"] %in% gStudio[,"ProbeID"]),]
gStudio <- gStudio[match(idatData[,"CodesBinData"], gStudio[,"ProbeID"]),]


###################################################
### code chunk number 8: figureComparingValues
###################################################
par(mfrow = c(1,2))
plot(idatData[, "MeanBinData"], gStudio[, "X4343238080_A.AVG_Signal"], 
     xlab = "illuminaio", ylab = "GenomeStudio")
identical(idatData[, "MeanBinData"], gStudio[, "X4343238080_A.AVG_Signal"])
hist(idatData[, "MeanBinData"]- gStudio[, "X4343238080_A.AVG_Signal"],
     breaks = 100, main = "", xlab = "Difference")


###################################################
### code chunk number 9: sessionInfo
###################################################
toLatex(sessionInfo())


