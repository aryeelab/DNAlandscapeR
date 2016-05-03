library(shiny)
source("global.R")

pageWithSidebar(
  headerPanel("DNAlandscapeR"),
  sidebarPanel(
    selectInput("tracks", label = h3("Tracks"), 
        choices = f.list, 
        selectize = TRUE, multiple= TRUE, selected = 0),
    textInput("chr", "Chromsome", value = "12", width = NULL, placeholder = NULL),
    textInput("start", "Start", value = "12660203", width = NULL, placeholder = NULL),
    textInput("stop", "Stop", value = "12777226", width = NULL, placeholder = NULL),
    textInput("Gene", "Gene", value = "AGO3", width = NULL, placeholder = NULL),
    checkboxInput("showgenes", "Show Gene Annotation", value = TRUE, width = NULL),
    actionButton("plot.region", "Plot Region"), 
    actionButton("plot.gene", "Plot Gene"),
    actionButton("zoom.in", "Zoom In"),
    actionButton("zoom.out", "Zoom Out"),
    downloadButton("down", "Download Plot"),
    actionButton("clear", "Clear")
  ),

  mainPanel(
      plotOutput("plot")
      )
)

