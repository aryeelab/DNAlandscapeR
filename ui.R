library(shiny)
source("global.R")

pageWithSidebar(
  headerPanel("DNA Regional Viewer"),
  sidebarPanel(
    selectInput("tracks", label = h3("Tracks"), 
        choices = f.list, 
        selectize = TRUE, multiple= TRUE, selected = 0),
    textInput("chr", "Chromsome", value = "1", width = NULL, placeholder = NULL),
    textInput("start", "Start", value = "36000000", width = NULL, placeholder = NULL),
    textInput("stop", "Stop", value = "36300000", width = NULL, placeholder = NULL),
    textInput("Gene", "Gene", value = "AGO3", width = NULL, placeholder = NULL),
    checkboxInput("showgenes", "Show Gene Annotation", value = TRUE, width = NULL),
    actionButton("plot.region", "Plot Region"), 
    actionButton("plot.gene", "Plot Gene"),
    downloadButton("down", "Download Plot"),
    actionButton("clear", "Clear")
  ),

  mainPanel(
      plotOutput("plot")
      )
)

