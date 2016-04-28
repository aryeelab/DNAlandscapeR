
library(shiny)
library(ggplot2)
library(GenomicRanges)
library(diffloop)
library(Sushi)

function(input, output) {
  val <- reactiveValues(region = NULL)
  
  rda<-paste(system.file('rda',package='diffloop'),'loops.small.rda',sep='/')
  load(rda)

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
    loopPlot(loops.small, val$region)
      }, height= 700)
 }
