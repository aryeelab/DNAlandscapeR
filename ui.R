# DNAlandscapeR UI # 
source("global.R")

#Custom input style
textInput3 <- function (inputId, label, value = "", ...){
    div(
    tags$label(label, `for` = inputId), 
    tags$input(id = inputId, type = "text", value = value, ...))
}

shinyUI(navbarPage(HTML("<img src='harvard-logo.png'/>"),
tabPanel("Welcome",
fluidPage(
    headerPanel(
        HTML("<h2><b><P ALIGN=Center>Welcome to the DNAlandscapeR Epigenomics Browser</b></h2>
             <h3><b><P ALIGN=Center>Aryee Lab</b></h3>")
        ),
    mainPanel(
        HTML("<h3><b>About</b></h3>"),
        tags$h5('We designed DNAlandscapeR to 
                provide a visualization of epigenetic data, particularly three-dimensional
                chromatin structures in a computationally efficient framework. In particular, 
                we utilize the R/Shiny environment to minimize the loaded data. If you use 
                DNAlandscapeR in your research, please cite <REFERENCE?>...'),
        width = 12
        )
)),
                   
tabPanel("Visualize", 
headerPanel(tags$h1(tags$b('DNA Landscape'))),
  sidebarLayout(
  sidebarPanel(
    uiOutput("trackoptions"), 
    tags$hr(),
    textInput3("chr", HTML("<h5><b>Chr&nbsp;&nbsp;&nbsp;&nbsp;</b></h5>"), value = default_chr),
    textInput3("start", HTML("<h5><b>Start&nbsp;</b></h5>"), value = default_start),
    textInput3("stop", HTML("<h5><b>Stop&nbsp;&nbsp;</b></h5>"), value = default_end),
    tags$br(),
    actionButton("plot.region", "Plot Region", style='padding:10px; font-size:80%'), 
    tags$hr(),
    actionButton("zoom.in", "Zoom In", style='padding:10px; font-size:80%'),
    actionButton("zoom.out", "Zoom Out", style='padding:10px; font-size:80%'),
    actionButton("clear", "Clear", style='padding:10px; font-size:80%'),
    tags$br(),
    tags$br(),
    actionButton("left.big", HTML("&nbsp; << &nbsp;"), style='padding:10px; font-size:80%'),
    actionButton("left.small", HTML("&nbsp; &nbsp; <  &nbsp; &nbsp;"), style='padding:10px; font-size:80%'),
    actionButton("right.small", HTML("&nbsp; &nbsp; > &nbsp; &nbsp;"), style='padding:10px; font-size:80%'),
    actionButton("right.big", HTML("&nbsp; >> &nbsp;"), style='padding:10px; font-size:80%'),
    tags$hr(),
    fileInput("skipRegions", HTML("<h4><b>Upload Regions .bed File</b></h4>"), multiple = FALSE, accept = NULL, width = NULL),
    actionButton("left.skip", "<<<", style='padding:10px; font-size:100%'),
    actionButton("right.skip", ">>>", style='padding:10px; font-size:100%'),
    tags$hr(),
    selectInput("showgenes", label=HTML("<h4><b>Genome Annotation</b></h4>"),
                choices = list("Gene Bodies" = 1, "Detailed Gene Annotation" = 2, "None" = 0), selected = 1)
  ), 

  mainPanel(
    #fluidRow(column(12, textOutput("regiontotal"), align = 'center')),
    plotOutput("plot", 
    dblclick = "plot_dblclick",
    brush = "plot_brush")
    ),
  fluid = TRUE
),
bsCollapse(id = "collapseAdvancedPlotOptions", open = "Panel1",
    bsCollapsePanel(title = HTML("<h4><b>Advanced options</b></h4>"), value = "Panel1",
    fluidRow(
       column(4, radioButtons("organism", HTML("<h4><b>Specify Organism</b></h4>"),
                    choices = list("Human" = 1, "Mouse" = 2), 
                    selected = 1),
              checkboxInput("log2BW", "log2 Transform .bigwigs", value = FALSE, width = NULL),
              checkboxInput("showSingleAnchors", "Show Single Anchors", value = FALSE, width = NULL)),
       column(4, textInput("Gene", HTML("<h4><b>Plot Gene Region </b></h4>"), value = "AGO3"),
                 actionButton("plot.gene", "Plot Gene", style='padding:10px; font-size:80%')),
       column(4, uiOutput("specifiedGenes"))
     ),
    tags$hr(),
    actionButton("refresh", "Refresh Plot"),
    downloadButton("down", "Download Plot"),    
    style = "default")
)),

tabPanel("Data Descriptions",
fluidPage(
    headerPanel(
        HTML("<h1><b>Preloaded Data Descriptions</b></h1>")
        ),
    mainPanel(
        tags$hr(),
        dataTableOutput('preloadedDataDescription'),
        tags$hr(), 
        width = 12
        )
)),

tabPanel("Upload",
pageWithSidebar(
    headerPanel(
        tags$h1(tags$b('Add local data'))
        ),
    sidebarPanel(
        radioButtons("organismUpload", HTML("<h5><b>Specify organism:</b></h5>"),
            choices = list("Human" = 1, "Mouse" = 2), 
            selected = 1),
        tags$hr(),
        fileInput("newTrack", h5(tags$b("Add file:")), multiple = FALSE, accept = NULL, width = NULL),
        tags$hr(),
        textInput3("newTrackName", h5(tags$b("Specify track name:"), value = "")),
        tags$hr(),
        selectInput("datType", label = h5(tags$b("Data type:")), 
                choices = uploadchoices, #defined in global.R
                            selected = 1, width = "90%"),
        actionButton("addFile", "Add", style='padding:10px; font-size:80%')
        ),
    mainPanel(
        HTML("<h4><b>Successfully loaded user files</b></h4>"),
        dataTableOutput('dt'),
        tags$hr()
        )
)),
tabPanel("Guide",
    includeMarkdown("www/DNAlandscapeR-help.Rmd")
),
theme = shinytheme("cosmo"),
footer = HTML('<P ALIGN=Center>DNAlandscapeR. &copy; <A HREF="mailto:caleblareau@g.harvard.edu">Caleb Lareau</A> & Martin Aryee'),
collapsible = TRUE, 
fluid = TRUE,
windowTitle = "DNAlandscapeR"
))
