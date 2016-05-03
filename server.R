library(shiny)
library(ggplot2)
library(GenomicRanges)
library(diffloop)
library(Sushi)
library(foreach)
library(rtracklayer)
library(tools)

source("helper.R")
source("global.R")

function(input, output, session) {
    val <- reactiveValues(region = NULL)
    observeEvent(input$plot.region, {
        val$region <- GRanges(seqnames=c(input$chr),ranges=IRanges(start=c(as.numeric(input$start)),end=c(as.numeric(input$stop))))
        })
  
    observeEvent(input$plot.gene, {
        rda<-paste(system.file('rda',package='diffloop'),'human.genes.rda',sep='/')
        load(rda)
        t <- human.genes[mcols(human.genes)$id == as.character(input$Gene) ]
        val$region <- padGRanges(t, pad = as.integer(width(t)/2))
        updateNumericInput(session, "chr", value = as.numeric(seqnames(val$region)))       
        updateNumericInput(session, "start", value = as.integer(start(ranges(range(val$region)))))  
        updateNumericInput(session, "stop", value = as.integer(end(ranges(range(val$region)))))       

    })  
  
    observeEvent(input$clear, {
        val$region <- NULL
    })


    output$plot <- renderPlot({
        if (is.null(val$region)) return()
        if (length(input$tracks) == 0) return()
        par(mfrow=c(length(input$tracks)+input$showgenes, 1), oma = c(0, 0, 0, 0), mar = c(3, 5, 1, 1))
        for(i in input$tracks){
            i <- as.integer(i)
            if (i < 1000){ #ChIA-PET from Data Source File
                oneSampleLoopPlot(paste("data/loops/", c.full[[i]], sep = ""),val$region)
            } else if (i < 2000) { # BigWig Read Count Track
                bw.trackplot(paste("data/tracks/", e.full[[i-1000]], sep = ""), val$region)
            } else if (i < 3000){ #DNA Methylation
                methyl.bedgraph.trackplot(paste("data/methylation/", m.full[[i-2000]],sep = ""), val$region)
            } else {return()}
        }
        if(input$showgenes) humanAnnotation(val$region)
     }, height = 700)
 }
