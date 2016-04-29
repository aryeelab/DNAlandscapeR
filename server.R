
library(shiny)
library(ggplot2)
library(GenomicRanges)
library(diffloop)
library(Sushi)
library(foreach)
library(gridExtra)

source("helper.R")
# Import data/annotate
rda<-paste(system.file('rda',package='diffloop'),'loops.small.rda',sep='/')
load(rda)
ctcf_j <- system.file('extdata','Jurkat_CTCF_chr1.narrowPeak',package='diffloop')
ctcf <- rmchr(padGRanges(bedToGRanges(ctcf_j), pad = 1000))
h3k27ac_j <- system.file('extdata','Jurkat_H3K27ac_chr1.narrowPeak',package='diffloop')
h3k27ac <- rmchr(padGRanges(bedToGRanges(h3k27ac_j), pad = 1000))
promoter <- padGRanges(getHumanTSS(c('1')), pad = 1000)
loops.small <- annotateLoops(loops.small, ctcf, h3k27ac, promoter)

function(input, output) {
    val <- reactiveValues(region = NULL)
    observeEvent(input$plot.region, {
        val$region <- GRanges(seqnames=c(input$chr),ranges=IRanges(start=c(as.numeric(input$start)),end=c(as.numeric(input$stop))))
        })
  
    observeEvent(input$plot.gene, {
        rda<-paste(system.file('rda',package='diffloop'),'human.genes.rda',sep='/')
        load(rda)
        t <- human.genes[mcols(human.genes)$id == as.character(input$Gene) ]
        val$region <- padGRanges(t, pad = as.integer(width(t)/2))
    })  
  
    observeEvent(input$clear, {
        val$region <- NULL
    })


    output$plot <- renderPlot({
        if (is.null(val$region)) return()
        if (length(input$tracks) == 0) return()
        par(mfrow=c(length(input$tracks)+1, 1), oma = c(0, 0, 0, 0), mar = c(3, 1, 1, 1))
        for(i in input$tracks){
            oneSampleLoopPlot(loops.small[,as.numeric(i)],val$region)
        }
        humanAnnotation(val$region)
     }, height = 700)
 }
