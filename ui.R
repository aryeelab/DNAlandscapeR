
source("global.R")

shinyUI(navbarPage("DNAlandscapeR",
tabPanel("Visualize", 
headerPanel('Visualize DNA Regions and Tracks'),
sidebarLayout(
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

  mainPanel(plotOutput("plot")),
  fluid = TRUE
)
),

tabPanel("Upload",
pageWithSidebar(
    headerPanel(
        'Add tracks from local files'
        ),
    sidebarPanel(
        tags$h4('Select additional file to be imported.'),
        shinyFiles::shinyFilesButton('file', 'File select', 'Please select a file', TRUE),
        tags$p(),
        tags$hr(),
        tags$p('The file selection button allows the user to select one or
               several files and get their absolute position communicated back
               to the shiny server. In this example the button has been set to
               single-file mode.'),
        radioButtons('datType', 'Data Type', c(Loops = 'loops', Methylation='methyl', Read.Depth="readdepth"), "readdepth"),
        radioButtons('fileformat', 'File Format', c(Bedgraph = 'bedgraph', BigWig='bigwig', rds="rds"), 'bigwig'),
        tags$hr()
        ),
    mainPanel(
        tags$h4('Current Specified File'),
        verbatimTextOutput('filein'),
        tags$hr()
        )
)
)))
