### R code from vignette source 'Sushi.Rnw'

###################################################
### code chunk number 1: Sushi.Rnw:25-27
###################################################
options(SweaveHooks=list(fig=function() par(mgp=c(3, .4, 0))))
options(continue=" ")


###################################################
### code chunk number 2: Sushi.Rnw:77-159
###################################################
writeLines("@article{encode_integrated_2012,
   title = {An integrated encyclopedia of {DNA} elements in the human genome},
   author = {{The ENCODE Project Consortium}},
   journal = {Nature},
   year = {2012},
   volume = {489},
   issn = {1476-4687},
   url = {http://www.ncbi.nlm.nih.gov/pubmed/22955616},
   doi = {10.1038/nature11247},
   number = {7414},
   month = sep,
   note = {{PMID:} 22955616},
   pages = {57--74},
 },
 @article{fiveC,
   title = {The long-range interaction landscape of gene promoters},
   volume = {489},
   doi = {10.1038/nature11279},
   number = {7414},
   journal = {Nature},
   author = {A Sanyal and BR Lajoie and G Jain and J Dekker},
   month = sep,
   year = {2012},
   note = {{PMID:} 22955621},
   pages = {109--113},
 },
 @article{chiapet,
   title = {Extensive promoter-centered chromatin interactions provide a topological basis for transcription regulation},
   volume = {148},
   doi = {10.1016/j.cell.2011.12.014},
   journal = {Cell},
   author = {G Li and X Ruan and RK Auerbach and KS Sandhu and M Zheng and P Wang and HM Poh and Y Goh and J Lim and J Zhang and HS Sim and SQ Peh and FH Mulawadi and CT Ong and YL Orlov and S Hong and Z Zhang and S Landt and D Raha and G Euskirchen and CL Wei and W Ge and H Wang and C Davis and KI Fisher-Aylor and A Mortazavi and M Gerstein and T Gingeras and B Wold and Y Sun and MJ Fullwood and E Cheung and E Liu and WK Sung and M Snyder and Y Ruan},
   month = jan,
   year = {2012},
   note = {{PMID:} 22265404},
   pages = {84--98},
 },
 @article{chipexo,
   title = {Comprehensive genome-wide protein-DNA interactions detected at single-nucleotide resolution},
   volume = {147},
   doi = {10.1016/j.cell.2011.11.013},
   journal = {Cell},
   author = {HS Rhee and BF Pugh},
   month = dec,
   year = {2011},
   note = {{PMID:} 22153082},
 },
 @article{dnaseI,
   title = {An expansive human regulatory lexicon encoded in transcription factor footprints},
   volume = {489},
   doi = {10.1038/nature11212},
   journal = {Nature},
   author = {S Neph and J Vierstra and AB Stergachis and AP Reynolds and E Haugen and B Vernot and RE Thurman and S John and R Sandstrom and AK Johnson and MT Maurano and R Humbert and E Rynes and H Wang and S Vong and K Lee and D Bates and M Diegel and V Roach and D Dunn and J Neri and A Schafer and RS Hansen and T  Kutyavin and E Giste and M Weaver and T Canfield and P Sabo and M Zhang and G Balasundaram and R Byron and MJ MacCoss and JM Akey and MA Bender and M Groudine and R Kaul and JA Stamatoyannopoulos},
   month = sep,
   year = {2012},
   pages = {83-90},
   note = {{PMID:} 22955618},
 },
 @article{GWAS,
   title = {Genetic variants in novel pathways influence blood pressure and cardiovascular disease risk},
   volume = {478},
   doi = {10.1038/nature10405},
   journal = {Nature},
   author = {{International Consortium for Blood Pressure}},
   month = sep,
   year = {2011},
   note = {{PMID:} 21909115},
 },
 @article{HiC,
  title = {Topological domains in mammalian genomes identified by analysis of chromatin interactions},
  volume = {485},
  doi = {10.1038/nature11082},
  journal = {Nature},
  author = {JR Dixon and S Selvaraj and F Yue and A Kim and Y Li and Y Shen and M Hu and JS Liu and B Ren},
  month = sep,
  year = {2012},
  note = {{PMID:} 22495300},
 },
 @online{biomart,
   author = {Biomart},
   url = {http://www.biomart.org/}
 }", con="Sushi.bib")


###################################################
### code chunk number 3: dataLoading
###################################################
library('Sushi')
Sushi_data = data(package = 'Sushi')
data(list = Sushi_data$results[,3]) 


###################################################
### code chunk number 4: dataLoading
###################################################
Sushi_data$results[,3]


###################################################
### code chunk number 5: Sushi.Rnw:202-203
###################################################
  head(Sushi_DNaseI.bedgraph)


###################################################
### code chunk number 6: Sushi.Rnw:211-215
###################################################
getOption("SweaveHooks")[["fig"]]()
chrom            = "chr11"
chromstart       = 1650000
chromend         = 2350000
plotBedgraph(Sushi_DNaseI.bedgraph,chrom,chromstart,chromend,colorbycol= SushiColors(5))


###################################################
### code chunk number 7: Sushi.Rnw:222-223 (eval = FALSE)
###################################################
## labelgenome(chrom,chromstart,chromend,n=4,scale="Mb")


###################################################
### code chunk number 8: Sushi.Rnw:226-231
###################################################
getOption("SweaveHooks")[["fig"]]()
chrom            = "chr11"
chromstart       = 1650000
chromend         = 2350000
plotBedgraph(Sushi_DNaseI.bedgraph,chrom,chromstart,chromend,colorbycol= SushiColors(5))
labelgenome(chrom,chromstart,chromend,n=4,scale="Mb")


###################################################
### code chunk number 9: Sushi.Rnw:240-242 (eval = FALSE)
###################################################
## mtext("Read Depth",side=2,line=1.75,cex=1,font=2)
## axis(side=2,las=2,tcl=.2)


###################################################
### code chunk number 10: Sushi.Rnw:245-252
###################################################
getOption("SweaveHooks")[["fig"]]()
chrom            = "chr11"
chromstart       = 1650000
chromend         = 2350000
plotBedgraph(Sushi_DNaseI.bedgraph,chrom,chromstart,chromend,colorbycol= SushiColors(5))
labelgenome(chrom,chromstart,chromend,n=4,scale="Mb")
mtext("Read Depth",side=2,line=1.75,cex=1,font=2)
axis(side=2,las=2,tcl=.2)


###################################################
### code chunk number 11: Sushi.Rnw:259-268
###################################################
getOption("SweaveHooks")[["fig"]]()
chrom            = "chr11"
chromstart       = 1955000
chromend         = 1960000
plotBedgraph(Sushi_ChIPSeq_CTCF.bedgraph,chrom,chromstart,chromend,
             transparency=.50,color=SushiColors(2)(2)[1])
plotBedgraph(Sushi_DNaseI.bedgraph,chrom,chromstart,chromend,
             transparency=.50,color=SushiColors(2)(2)[2],overlay=TRUE,
             rescaleoverlay=TRUE)
labelgenome(chrom,chromstart,chromend,n=3,scale="Kb")


###################################################
### code chunk number 12: Sushi.Rnw:276-280 (eval = FALSE)
###################################################
## 
## legend("topright",inset=0.025,legend=c("DNaseI","ChIP-seq (CTCF)"),
##        fill=opaque(SushiColors(2)(2)),border=SushiColors(2)(2),text.font=2,
##        cex=1.0)


###################################################
### code chunk number 13: Sushi.Rnw:282-295
###################################################
getOption("SweaveHooks")[["fig"]]()
chrom            = "chr11"
chromstart       = 1955000
chromend         = 1960000
plotBedgraph(Sushi_ChIPSeq_CTCF.bedgraph,chrom,chromstart,chromend,
             transparency=.50,color=SushiColors(2)(2)[1])
plotBedgraph(Sushi_DNaseI.bedgraph,chrom,chromstart,chromend,
             transparency=.50,color=SushiColors(2)(2)[2],overlay=TRUE,
             rescaleoverlay=TRUE)
labelgenome(chrom,chromstart,chromend,n=3,scale="Kb")

legend("topright",inset=0.025,legend=c("DNaseI","ChIP-seq (CTCF)"),
       fill=opaque(SushiColors(2)(2)),border=SushiColors(2)(2),text.font=2,
       cex=1.0)


###################################################
### code chunk number 14: Sushi.Rnw:303-304 (eval = FALSE)
###################################################
## par(mfrow=c(2,1),mar=c(1,4,1,1))


###################################################
### code chunk number 15: Sushi.Rnw:311-318 (eval = FALSE)
###################################################
## plotBedgraph(Sushi_ChIPSeq_CTCF.bedgraph,chrom,chromstart,chromend,transparency=.50,
##              color=SushiColors(2)(2)[1])
## axis(side=2,las=2,tcl=.2)
## mtext("Read Depth",side=2,line=1.75,cex=1,font=2)
## legend("topright",inset=0.025,legend=c("DNaseI","ChIP-seq (CTCF)"),
##        fill=opaque(SushiColors(2)(2)),border=SushiColors(2)(2),text.font=2,
##        cex=1.0)


###################################################
### code chunk number 16: Sushi.Rnw:324-340
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mfrow=c(2,1),mar=c(1,4,1,1))
# set the genomic regions
chrom            = "chr11"
chromstart       = 1955000
chromend         = 1960000

# plot chip-seq data
plotBedgraph(Sushi_ChIPSeq_CTCF.bedgraph,chrom,chromstart,chromend,transparency=.50)

# add y-axis
axis(side=2,las=2,tcl=.2)
mtext("Read Depth",side=2,line=1.75,cex=1,font=2)

# add legend
legend("topright",inset=0.025,legend=c("DNaseI","ChIP-seq (CTCF)"),
       fill=opaque(SushiColors(2)(2)),border=SushiColors(2)(2),text.font=2)


###################################################
### code chunk number 17: Sushi.Rnw:347-353 (eval = FALSE)
###################################################
## plotBedgraph(Sushi_DNaseI.bedgraph, chrom, chromstart, chromend,
##              transparency=.50, flip=TRUE, color=SushiColors(2)(2)[2])
## labelgenome(chrom,chromstart,chromend,side=3,n=3,scale="Kb")
## axis(side=2,las=2,tcl=.2,at=pretty(par("yaxp")[c(1,2)]),
##              labels=-1*pretty(par("yaxp")[c(1,2)]))
## mtext("Read Depth",side=2,line=1.75,cex=1,font=2)


###################################################
### code chunk number 18: Sushi.Rnw:359-388
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mfrow=c(2,1),mar=c(1,4,1,1))
# set the genomic regions
chrom            = "chr11"
chromstart       = 1955000
chromend         = 1960000

# plot chip-seq data
plotBedgraph(Sushi_ChIPSeq_CTCF.bedgraph ,chrom, chromstart, chromend,
             transparency=.50, color=SushiColors(2)(2)[1])

# add y-axis
axis(side=2,las=2,tcl=.2)
mtext("Read Depth",side=2,line=1.75,cex=1,font=2)

# add legend
legend("topright",inset=0.025,legend=c("DNaseI","ChIP-seq (CTCF)"),
       fill=opaque(SushiColors(2)(2)),border=SushiColors(2)(2),text.font=2)

# plot dnaseI data
plotBedgraph(Sushi_DNaseI.bedgraph,chrom,chromstart,chromend,transparency=.50,flip=TRUE,
             color=SushiColors(2)(2)[2])

# add y-axis
ylabs = axis(side=2,las=2,tcl=.2,at=pretty(par("yaxp")[c(1,2)]),
             labels=-1*pretty(par("yaxp")[c(1,2)]))
mtext("Read Depth",side=2,line=1.75,cex=1,font=2)

# add the genome labels
labelgenome(chrom,chromstart,chromend,side=3,n=3,scale="Kb")


###################################################
### code chunk number 19: Sushi.Rnw:397-398
###################################################
  Sushi_HiC.matrix[100:105,100:105]


###################################################
### code chunk number 20: Sushi.Rnw:405-415
###################################################
getOption("SweaveHooks")[["fig"]]()
chrom            = "chr11"
chromstart       = 500000
chromend         = 5050000
phic = plotHic(Sushi_HiC.matrix, chrom,chromstart, chromend, max_y = 20,
               zrange=c(0,28), palette=SushiColors(7))
addlegend(phic[[1]], palette=phic[[2]], title="score", side="right",
          bottominset=0.4, topinset=0, xoffset=-.035, labelside="left",
          width=0.025, title.offset=0.035)
labelgenome(chrom, chromstart, chromend, n=4, scale="Mb",
            edgeblankfraction=0.20)


###################################################
### code chunk number 21: Sushi.Rnw:427-438
###################################################
getOption("SweaveHooks")[["fig"]]()
chrom            = "chr11"
chromstart       = 500000
chromend         = 5050000

phic = plotHic(Sushi_HiC.matrix,chrom,chromstart,chromend,max_y = 20,
               zrange=c(0,28),flip=TRUE,palette=topo.colors)

addlegend(phic[[1]],palette=phic[[2]],title="score",side="left",bottominset=0.1,
          topinset=0.5,xoffset=-.035,labelside="right",width=0.025,title.offset=0.035)

labelgenome(chrom,chromstart,chromend,side=3,n=4,scale="Mb",edgeblankfraction=0.20)


###################################################
### code chunk number 22: Sushi.Rnw:448-449
###################################################
  head(Sushi_5C.bedpe)


###################################################
### code chunk number 23: Sushi.Rnw:456-468
###################################################
getOption("SweaveHooks")[["fig"]]()
chrom            = "chr11"
chromstart       = 1650000
chromend         = 2350000
pbpe = plotBedpe(Sushi_5C.bedpe,chrom,chromstart,chromend,
                 heights = Sushi_5C.bedpe$score,plottype="loops",
                 colorby=Sushi_5C.bedpe$samplenumber,
                 colorbycol=SushiColors(3))
labelgenome(chrom, chromstart,chromend,n=3,scale="Mb")
legend("topright",inset =0.01,legend=c("K562","HeLa","GM12878"),
       col=SushiColors(3)(3),pch=19,bty='n',text.font=2)
axis(side=2,las=2,tcl=.2)
mtext("Z-score",side=2,line=1.75,cex=.75,font=2)


###################################################
### code chunk number 24: Sushi.Rnw:476-485
###################################################
getOption("SweaveHooks")[["fig"]]()
chrom            = "chr11"
chromstart       = 1650000
chromend         = 2350000
pbpe = plotBedpe(Sushi_5C.bedpe,chrom,chromstart,chromend,flip=TRUE,
                 plottype="lines",colorby=Sushi_5C.bedpe$score,
                 colorbycol=SushiColors(5))
labelgenome(chrom, chromstart,chromend,side=3,n=3,scale="Mb")
addlegend(pbpe[[1]],palette=pbpe[[2]],title="Z-score",side="right",bottominset=0.05,
          topinset=0.05,xoffset=-.035,labelside="right",width=0.025,title.offset=0.045)


###################################################
### code chunk number 25: Sushi.Rnw:495-496
###################################################
  head(Sushi_ChIPSeq_pol2.bed)


###################################################
### code chunk number 26: Sushi.Rnw:502-514
###################################################
getOption("SweaveHooks")[["fig"]]()
chrom            = "chr11"
chromstart       = 2281200
chromend         = 2282200

plotBed(beddata    = Sushi_ChIPSeq_pol2.bed,chrom = chrom,chromstart = chromstart,
        chromend =chromend,colorby    = Sushi_ChIPSeq_pol2.bed$strand,
        colorbycol = SushiColors(2),row  = "auto",wiggle=0.001)

labelgenome(chrom,chromstart,chromend,n=2,scale="Kb")

legend("topright",inset=0,legend=c("reverse","forward"),fill=SushiColors(2)(2),
       border=SushiColors(2)(2),text.font=2,cex=0.75)


###################################################
### code chunk number 27: Sushi.Rnw:521-533
###################################################
getOption("SweaveHooks")[["fig"]]()
chrom            = "chr11"
chromstart       = 2281200
chromend         = 2282200

plotBed(beddata    = Sushi_ChIPSeq_pol2.bed,chrom = chrom,chromstart = chromstart,
        chromend =chromend,colorby    = Sushi_ChIPSeq_pol2.bed$strand,
        colorbycol = SushiColors(2),row  = "auto",wiggle=0.001,splitstrand=TRUE)

labelgenome(chrom,chromstart,chromend,n=2,scale="Kb")

legend("topright",inset=0,legend=c("reverse","forward"),fill=SushiColors(2)(2),
       border=SushiColors(2)(2),text.font=2,cex=0.75)


###################################################
### code chunk number 28: Sushi.Rnw:540-543
###################################################
Sushi_ChIPSeq_severalfactors.bed$color = 
        maptocolors(Sushi_ChIPSeq_severalfactors.bed$row,
        col=SushiColors(6))


###################################################
### code chunk number 29: Sushi.Rnw:550-568
###################################################
getOption("SweaveHooks")[["fig"]]()
chrom            = "chr15"
chromstart      = 72800000
chromend         = 73100000

Sushi_ChIPSeq_severalfactors.bed$color = 
        maptocolors(Sushi_ChIPSeq_severalfactors.bed$row,
        col=SushiColors(6))

plotBed(beddata    = Sushi_ChIPSeq_severalfactors.bed,chrom = chrom,
        chromstart = chromstart,chromend =chromend,
        rownumber  = Sushi_ChIPSeq_severalfactors.bed$row, type = "circles",
        color=Sushi_ChIPSeq_severalfactors.bed$color,row="given",
        plotbg="grey95",rowlabels=unique(Sushi_ChIPSeq_severalfactors.bed$name),
        rowlabelcol=unique(Sushi_ChIPSeq_severalfactors.bed$color),rowlabelcex=0.75)

labelgenome(chrom,chromstart,chromend,n=3,scale="Mb")

mtext("ChIP-seq",side=3, adj=-0.065,line=0.5,font=2)


###################################################
### code chunk number 30: Sushi.Rnw:574-588
###################################################
chrom            = "chr15"
chromstart      = 72800000
chromend         = 73100000

plotBed(beddata    = Sushi_ChIPSeq_severalfactors.bed,chrom = chrom,
        chromstart = chromstart,chromend =chromend,
        rownumber  = Sushi_ChIPSeq_severalfactors.bed$row, type = "circles",
        color=Sushi_ChIPSeq_severalfactors.bed$color,row="given",
        plotbg="grey95",rowlabels=unique(Sushi_ChIPSeq_severalfactors.bed$name),
        rowlabelcol=unique(Sushi_ChIPSeq_severalfactors.bed$color),rowlabelcex=0.75)

labelgenome(chrom,chromstart,chromend,n=3,scale="Mb")

mtext("ChIP-seq",side=3, adj=-0.065,line=0.5,font=2)


###################################################
### code chunk number 31: Sushi.Rnw:595-605
###################################################
plotBed(beddata    = Sushi_ChIPSeq_severalfactors.bed,chrom = chrom,
        chromstart = chromstart,chromend =chromend,
        rownumber  = Sushi_ChIPSeq_severalfactors.bed$row, type = "region",
        color=Sushi_ChIPSeq_severalfactors.bed$color,row="given",
        plotbg="grey95",rowlabels=unique(Sushi_ChIPSeq_severalfactors.bed$name),
        rowlabelcol=unique(Sushi_ChIPSeq_severalfactors.bed$color),rowlabelcex=0.75)

labelgenome(chrom,chromstart,chromend,n=3,scale="Mb")

mtext("ChIP-seq",side=3, adj=-0.065,line=0.5,font=2)


###################################################
### code chunk number 32: Sushi.Rnw:610-628
###################################################
getOption("SweaveHooks")[["fig"]]()
chrom            = "chr15"
chromstart       = 72800000
chromend         = 73100000

Sushi_ChIPSeq_severalfactors.bed$color = 
        maptocolors(Sushi_ChIPSeq_severalfactors.bed$row,
        col=SushiColors(6))

plotBed(beddata    = Sushi_ChIPSeq_severalfactors.bed,chrom = chrom,
        chromstart = chromstart,chromend =chromend,
        rownumber  = Sushi_ChIPSeq_severalfactors.bed$row, type = "region",
        color=Sushi_ChIPSeq_severalfactors.bed$color,row="given",
        plotbg="grey95",rowlabels=unique(Sushi_ChIPSeq_severalfactors.bed$name),
        rowlabelcol=unique(Sushi_ChIPSeq_severalfactors.bed$color),rowlabelcex=0.75)

labelgenome(chrom,chromstart,chromend,n=3,scale="Mb")

mtext("ChIP-seq",side=3, adj=-0.065,line=0.5,font=2)


###################################################
### code chunk number 33: Sushi.Rnw:635-648
###################################################
chrom            = "chr15"
chromstart       = 60000000
chromend         = 80000000
chrom_biomart    = gsub("chr","",chrom)

mart=useMart(host='may2009.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', 
             dataset='hsapiens_gene_ensembl')

geneinfobed = getBM(attributes = c("chromosome_name","start_position","end_position"),
                    filters= c("chromosome_name","start","end"),
                    values=list(chrom_biomart,chromstart,chromend),mart=mart)

geneinfobed[,1] = paste("chr",geneinfobed[,1],sep="")


###################################################
### code chunk number 34: Sushi.Rnw:655-656
###################################################
head (geneinfobed)


###################################################
### code chunk number 35: Sushi.Rnw:664-671
###################################################
plotBed(beddata = geneinfobed[!duplicated(geneinfobed),],chrom = chrom,
        chromstart = chromstart,chromend =chromend,row='supplied',
        palettes = list(SushiColors(7)), type = "density")

labelgenome(chrom, chromstart, chromend,  n=4,scale="Mb",edgeblankfraction=0.10)

mtext("Gene Density",side=3, adj=0,line=0.20,font=2)


###################################################
### code chunk number 36: Sushi.Rnw:676-698
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mar=c(3,1,3,1))
chrom            = "chr15"
chromstart       = 60000000
chromend         = 80000000
chrom_biomart    = gsub("chr","",chrom)

mart=useMart(host='may2009.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', 
             dataset='hsapiens_gene_ensembl')

geneinfobed = getBM(attributes = c("chromosome_name","start_position","end_position"),
                    filters= c("chromosome_name","start","end"),
                    values=list(chrom_biomart,chromstart,chromend),mart=mart)

geneinfobed[,1] = paste("chr",geneinfobed[,1],sep="")

plotBed(beddata = geneinfobed[!duplicated(geneinfobed),],chrom = chrom,
        chromstart = chromstart,chromend =chromend,row='supplied',
        palettes = list(SushiColors(7)), type = "density")

labelgenome(chrom, chromstart, chromend,  n=4,scale="Mb",edgeblankfraction=0.10)

mtext("Gene Density",side=3, adj=0,line=0.20,font=2)


###################################################
### code chunk number 37: Sushi.Rnw:711-712
###################################################
  head(Sushi_GWAS.bed)


###################################################
### code chunk number 38: Sushi.Rnw:717-718
###################################################
  head(Sushi_hg18_genome)


###################################################
### code chunk number 39: Sushi.Rnw:724-731
###################################################
plotManhattan(bedfile=Sushi_GWAS.bed,pvalues=Sushi_GWAS.bed[,5],
                col=SushiColors(6),genome=Sushi_hg18_genome,cex=0.75)
labelgenome(genome=Sushi_hg18_genome,n=4,scale="Mb",
                edgeblankfraction=0.20,cex.axis=.5)
axis(side=2,las=2,tcl=.2)
mtext("log10(P)",side=2,line=1.75,cex=1,font=2)
mtext("chromosome",side=1,line=1.75,cex=1,font=2)


###################################################
### code chunk number 40: GWAS1
###################################################
png('GWAS1.png',height=600,width=800 )
plotManhattan(bedfile=Sushi_GWAS.bed,pvalues=Sushi_GWAS.bed[,5],
                col=SushiColors(6),genome=Sushi_hg18_genome,cex=0.75)
labelgenome(genome=Sushi_hg18_genome,n=4,scale="Mb",
                edgeblankfraction=0.20,cex.axis=.5)
axis(side=2,las=2,tcl=.2)
mtext("log10(P)",side=2,line=1.75,cex=1,font=2)
mtext("chromosome",side=1,line=1.75,cex=1,font=2)
dev.off() 


###################################################
### code chunk number 41: Sushi.Rnw:755-756
###################################################
head(Sushi_genes.bed)


###################################################
### code chunk number 42: Sushi.Rnw:761-772
###################################################
getOption("SweaveHooks")[["fig"]]()
chrom            = "chr15"
chromstart       = 72998000
chromend         = 73020000

pg = plotGenes(Sushi_genes.bed,chrom,chromstart,chromend ,
               types=Sushi_genes.bed$type,maxrows=1,bheight=0.2,
               plotgenetype="arrow",bentline=FALSE,
          labeloffset=.4,fontsize=1.2,arrowlength = 0.025,
               labeltext=TRUE)

labelgenome( chrom, chromstart,chromend,n=3,scale="Mb")


###################################################
### code chunk number 43: Sushi.Rnw:778-779
###################################################
Sushi_transcripts.bed[1:20,]


###################################################
### code chunk number 44: Sushi.Rnw:786-802
###################################################
getOption("SweaveHooks")[["fig"]]()
chrom            = "chr15"
chromstart       = 72965000
chromend         = 72990000

pg = plotGenes(Sushi_transcripts.bed,chrom,chromstart,chromend ,
               types = Sushi_transcripts.bed$type,
               colorby=log10(Sushi_transcripts.bed$score+0.001),
               colorbycol= SushiColors(5),colorbyrange=c(0,1.0),
               labeltext=TRUE,maxrows=50,height=0.4,plotgenetype="box")


labelgenome( chrom, chromstart,chromend,n=3,scale="Mb")

addlegend(pg[[1]],palette=pg[[2]],title="log10(FPKM)",side="right",
          bottominset=0.4,topinset=0,xoffset=-.035,labelside="left",
          width=0.025,title.offset=0.055)


###################################################
### code chunk number 45: Sushi.Rnw:816-818 (eval = FALSE)
###################################################
## layout(matrix(c(1,1,2,3),2, 2, byrow = TRUE))
## par(mar=c(3,4,1,1))


###################################################
### code chunk number 46: Sushi.Rnw:823-835 (eval = FALSE)
###################################################
## chrom            = "chr11"
## chromstart       = 1900000
## chromend         = 2350000
## 
## plotBedgraph(Sushi_DNaseI.bedgraph,chrom,chromstart=chromstart,
##              chromend=chromend,colorbycol= SushiColors(5))
## 
## labelgenome(chrom,chromstart=chromstart,chromend=chromend,n=4,
##             scale="Mb")
## 
## mtext("Read Depth",side=2,line=1.75,cex=1,font=2)
## axis(side=2,las=2,tcl=.2)


###################################################
### code chunk number 47: Sushi.Rnw:840-851
###################################################
getOption("SweaveHooks")[["fig"]]()
layout(matrix(c(1,1,2,3),2, 2, byrow = TRUE))
par(mar=c(3,4,1,1))
chrom            = "chr11"
chromstart       = 1900000
chromend         = 2350000
plotBedgraph(Sushi_DNaseI.bedgraph,chrom,chromstart=chromstart,
             chromend=chromend,colorbycol= SushiColors(5))
labelgenome(chrom,chromstart=chromstart,chromend=chromend,n=4,
            scale="Mb")
mtext("Read Depth",side=2,line=1.75,cex=1,font=2)
axis(side=2,las=2,tcl=.2)


###################################################
### code chunk number 48: Sushi.Rnw:857-864 (eval = FALSE)
###################################################
## zoomregion1      = c(1955000,1960000)
## zoomregion2      = c(2279000,2284000)
## 
## zoomsregion(zoomregion1,extend=c(0.01,0.13),wideextend=0.05,
##             offsets=c(0,0.580))
## zoomsregion(zoomregion2,extend=c(0.01,0.13),wideextend=0.05,
##             offsets=c(0.580,0))


###################################################
### code chunk number 49: Sushi.Rnw:868-887
###################################################
getOption("SweaveHooks")[["fig"]]()
layout(matrix(c(1,1,2,3),2, 2, byrow = TRUE))
par(mar=c(3,4,1,1))

chrom            = "chr11"
chromstart       = 1900000
chromend         = 2350000

plotBedgraph(Sushi_DNaseI.bedgraph,chrom,chromstart=chromstart,
             chromend=chromend,colorbycol= SushiColors(5))

labelgenome(chrom,chromstart=chromstart,chromend=chromend,n=4,scale="Mb")

mtext("Read Depth",side=2,line=1.75,cex=1,font=2)
axis(side=2,las=2,tcl=.2)

zoomregion1      = c(1955000,1960000)
zoomregion2      = c(2279000,2284000)
zoomsregion(zoomregion1,extend=c(0.01,0.13),wideextend=0.05,offsets=c(0,0.580))
zoomsregion(zoomregion2,extend=c(0.01,0.13),wideextend=0.05,offsets=c(0.580,0))


###################################################
### code chunk number 50: Sushi.Rnw:893-910 (eval = FALSE)
###################################################
## plotBedgraph(Sushi_DNaseI.bedgraph,chrom,chromstart=zoomregion1[1],
##              chromend=zoomregion1[2],colorbycol= SushiColors(5))
## 
## labelgenome(chrom,chromstart=zoomregion1[1],chromend=zoomregion1[2],
##             n=4,scale="Kb",edgeblankfraction=0.2,cex.axis=.75)
## zoombox()
## mtext("Read Depth",side=2,line=1.75,cex=1,font=2)
## axis(side=2,las=2,tcl=.2)
## 
## plotBedgraph(Sushi_DNaseI.bedgraph,chrom,chromstart=zoomregion2[1],
##              chromend=zoomregion2[2],colorbycol= SushiColors(5))
## 
## labelgenome(chrom,chromstart=zoomregion2[1],chromend=zoomregion2[2],
##             n=4,scale="Kb",edgeblankfraction=0.2,cex.axis=.75)
## zoombox()
## mtext("Read Depth",side=2,line=1.75,cex=1,font=2)
## axis(side=2,las=2,tcl=.2)


###################################################
### code chunk number 51: Sushi.Rnw:915-958
###################################################
getOption("SweaveHooks")[["fig"]]()
layout(matrix(c(1,1,2,3),2, 2, byrow = TRUE))
par(mar=c(3,4,1,1))

chrom            = "chr11"
chromstart       = 1900000
chromend         = 2350000

plotBedgraph(Sushi_DNaseI.bedgraph,chrom,chromstart=chromstart,
             chromend=chromend,colorbycol= SushiColors(5))

labelgenome(chrom,chromstart=chromstart,chromend=chromend,n=4,
            scale="Mb")

mtext("Read Depth",side=2,line=1.75,cex=1,font=2)
axis(side=2,las=2,tcl=.2)

zoomregion1      = c(1955000,1960000)
zoomregion2      = c(2279000,2284000)
zoomsregion(zoomregion1,extend=c(0.01,0.13),wideextend=0.05,
            offsets=c(0,0.580))
zoomsregion(zoomregion2,extend=c(0.01,0.13),wideextend=0.05,
            offsets=c(0.580,0))

plotBedgraph(Sushi_DNaseI.bedgraph,chrom,chromstart=zoomregion1[1],
             chromend=zoomregion1[2],colorbycol= SushiColors(5))

labelgenome(chrom,chromstart=zoomregion1[1],chromend=zoomregion1[2],
            n=4,scale="Kb",edgeblankfraction=0.2,cex.axis=.75)
zoombox()

mtext("Read Depth",side=2,line=1.75,cex=1,font=2)
axis(side=2,las=2,tcl=.2)

plotBedgraph(Sushi_DNaseI.bedgraph,chrom,chromstart=zoomregion2[1],
             chromend=zoomregion2[2],colorbycol= SushiColors(5))

labelgenome(chrom,chromstart=zoomregion2[1],chromend=zoomregion2[2],
            n=4,scale="Kb",edgeblankfraction=0.2,cex.axis=.75)
zoombox()

mtext("Read Depth",side=2,line=1.75,cex=1,font=2)
axis(side=2,las=2,tcl=.2)



###################################################
### code chunk number 52: Sushi.Rnw:975-976
###################################################
SushiColors(palette='list')


###################################################
### code chunk number 53: Sushi.Rnw:981-996
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(1,xlab='',xaxt='n',ylab='',yaxt='n',xlim=c(0.5,7.5),
     ylim=c(2,7.5),type='n')
for (i in (2:7))
{
  for (j in (1:i))
  {
    rect(j-.5,i,j+.5,i+.5,col=SushiColors(i)(i)[j])
  }
}

axis(side=2,at=(2:7),labels=(2:7),las=2)
axis(side=1,at=(1:7),labels=(1:7))
mtext("SushiColors",side=3,font=2, line=1, cex=1.5)
mtext("colors",side=1,font=2, line=2)
mtext("palette",side=2,font=2, line=2)


###################################################
### code chunk number 54: Sushi.Rnw:1003-1017
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(1,xlab='',xaxt='n',ylab='',yaxt='n',bty='n',type='n',
     xlim=c(-.15,1.05),ylim=c(-1,2))
for (i in seq(0,1,by=0.1))
{
    rect(i-.05,-1,i+.05,1,col=opaque("red",transparency=i))
    rect(i-.05,0,i+.05,2,col=opaque("blue",transparency=1-i))
}
axis(side=1,at=seq(0,1,by=0.1),labels=seq(0,1,by=0.1))
mtext("red transparency",side=1,font=2, line=2)
axis(side=3,at=seq(0,1,by=0.1),labels=seq(1,0,by=-0.1))
mtext("blue transparency",side=3,font=2, line=2)
text(-0.075,1.5,labels="blue",font=2,adj=1)
text(-0.075,0.5,labels="overlap",font=2,adj=1)
text(-0.075,-.5,labels="red",font=2,adj=1)


###################################################
### code chunk number 55: Sushi.Rnw:1024-1032
###################################################
getOption("SweaveHooks")[["fig"]]()
set.seed(3)
values = rnorm((1:10))
colorpalette = SushiColors(5)
plot(x=(1:10),y=values,col=maptocolors(values,colorpalette),
     pch=19,cex=4,xlab="data points",yaxt='n',ylim=range(values)*1.2)
addlegend(range(values),title="key",palette=colorpalette,
          side='left',xoffset = -0.125,width=0.03,bottominset = 0.5, topinset = 0.025)
axis(side=2,las=2)


###################################################
### code chunk number 56: Sushi.Rnw:1043-1053
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mar=c(8,3,3,1),mgp=c(3, .3, 0))

plotBedgraph(Sushi_DNaseI.bedgraph,chrom="chr11",chromstart=1650000,
             chromend=2350000,colorbycol=SushiColors(7))
labelgenome(chrom="chr11",chromstart=1650000,chromend=2350000,
            side=1,n=4,scale="Mb",line=.25)
labelgenome(chrom="chr11",chromstart=1650000,chromend=2350000,
            side=1,n=3,scale="Kb",line=2)
labelgenome(chrom="chr11",chromstart=1650000,chromend=2350000,
            side=1,n=1,scale="bp",line=4)


###################################################
### code chunk number 57: Sushi.Rnw:1060-1068
###################################################
plotManhattan(bedfile=Sushi_GWAS.bed,pvalues=Sushi_GWAS.bed[,5],
              col=SushiColors(6),genome=Sushi_hg18_genome,
              cex=0.75,space=0.05)
labelgenome(genome=Sushi_hg18_genome,n=4,scale="Mb",
            edgeblankfraction=0.20,cex.axis=.5,space=0.05)
axis(side=2,las=2,tcl=.2)
mtext("log10(P)",side=2,line=1.75,cex=1,font=2)
mtext("chromosome",side=1,line=1.75,cex=1,font=2)


###################################################
### code chunk number 58: GWAS2
###################################################
png('GWAS2.png',height=600,width=900 )
plotManhattan(bedfile=Sushi_GWAS.bed,pvalues=Sushi_GWAS.bed[,5],
              col=SushiColors(6),genome=Sushi_hg18_genome,
              cex=0.75,space=0.05)
labelgenome(genome=Sushi_hg18_genome,n=4,scale="Mb",
            edgeblankfraction=0.20,cex.axis=.5,space=0.05)
axis(side=2,las=2,tcl=.2)
mtext("log10(P)",side=2,line=1.75,cex=1,font=2)
mtext("chromosome",side=1,line=1.75,cex=1,font=2)
dev.off()


###################################################
### code chunk number 59: Sushi.Rnw:1092-1093
###################################################
labelplot("A) ","Manhattan Plot")


###################################################
### code chunk number 60: Sushi.Rnw:1097-1106
###################################################
plotManhattan(bedfile=Sushi_GWAS.bed,pvalues=Sushi_GWAS.bed[,5],
              col=SushiColors(6),genome=Sushi_hg18_genome,
              cex=0.75,space=0.05)
labelgenome(genome=Sushi_hg18_genome,n=4,scale="Mb"
            ,edgeblankfraction=0.20,cex.axis=.5,space=0.05)
axis(side=2,las=2,tcl=.2)
mtext("log10(P)",side=2,line=1.75,cex=1,font=2)
mtext("chromosome",side=1,line=1.75,cex=1,font=2)
labelplot("A) ","Manhattan Plot")


###################################################
### code chunk number 61: GWAS3
###################################################
png('GWAS3.png',height=600,width=900 )
plotManhattan(bedfile=Sushi_GWAS.bed,pvalues=Sushi_GWAS.bed[,5],
              col=SushiColors(6),genome=Sushi_hg18_genome,
              cex=0.75,space=0.05)
labelgenome(genome=Sushi_hg18_genome,n=4,scale="Mb"
            ,edgeblankfraction=0.20,cex.axis=.5,space=0.05)
axis(side=2,las=2,tcl=.2)
mtext("log10(P)",side=2,line=1.75,cex=1,font=2)
mtext("chromosome",side=1,line=1.75,cex=1,font=2)
labelplot("A) ","Manhattan Plot")
dev.off()


###################################################
### code chunk number 62: Sushi.Rnw:1126-1130
###################################################
if (file.exists("Rplots.pdf"))
{
  file.remove("Rplots.pdf")
}


###################################################
### code chunk number 63: Sushi.Rnw:1160-1161 (eval = FALSE)
###################################################
## read.table(file="reads.bed",sep="\t")


