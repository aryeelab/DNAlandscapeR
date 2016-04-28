
library(shiny)
library(ggplot2)
library(GenomicRanges)
library(diffloop)
library(Sushi)

pageWithSidebar(
  headerPanel("DNA Regional Viewer"),
  sidebarPanel(
    textInput("chr", "Chromsome", value = "1", width = NULL, placeholder = NULL),
    textInput("start", "Start", value = "36000000", width = NULL, placeholder = NULL),
    textInput("stop", "Stop", value = "36100000", width = NULL, placeholder = NULL),
    textInput("Gene", "Gene", value = "MYC", width = NULL, placeholder = NULL),
    actionButton("plot.region", "Plot Region"), 
    actionButton("plot.gene", "Plot Gene")
  ),

  mainPanel(plotOutput("plot"))
)

