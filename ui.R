
source("global.R")

shinyUI(navbarPage("DNAlandscapeR",
tabPanel("Visualize", 
headerPanel('Visualize DNA Regions and Tracks'),
sidebarLayout(
  sidebarPanel(
    uiOutput("trackoptions"), 
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

  mainPanel(plotOutput("plot")),
  fluid = TRUE
)
),

tabPanel("Upload",
pageWithSidebar(
    headerPanel(
        'Add local tracks'
        ),
    sidebarPanel(
        tags$h4('Choose File'),
        verbatimTextOutput('filename'),
        shinyFilesButton('file', 'File select', 'Please select a file', TRUE),
        tags$hr(),
        tags$h4('Specify input data formats'),
        radioButtons('datType', 'Data Type', c(Loops = 'Loops', Methylation='Methyl', Read.Depth="Read.Depth"), "Loops"),
        radioButtons('fileformat', 'File Format', c(Bedgraph = 'Bedgraph', BigWig='BigWig', rds="rds"), 'rds'),
        tags$hr(),
        actionButton("addFile", "Add")
        ),
    mainPanel(
        tags$h4('Data Frame of successfuly loaded user files'),
        dataTableOutput('dt'),
        tags$hr()
        )
)
)))
