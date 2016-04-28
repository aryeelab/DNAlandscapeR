
library(shiny)
library(ggplot2)
library(GenomicRanges)
library(diffloop)
library(Sushi)

function(input, output) {

  rda<-paste(system.file('rda',package='diffloop'),'loops.small.rda',sep='/')
  load(rda)
    
  region <- eventReactive(input$plot.region, {
     GRanges(seqnames=c(input$chr),ranges=IRanges(start=c(as.numeric(input$start)),end=c(as.numeric(input$stop))))
  })
  
  output$plot <- renderPlot({loopPlot(loops.small, region())}, height= 700)
    
 }
