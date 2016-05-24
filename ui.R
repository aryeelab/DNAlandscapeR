
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
        tags$h3('About:'),
        tags$h6('We designed DNAlandscapeR to 
                provide a visualization of epigenetic data, particularly three-dimensional
                chromatin structures in a computationally efficient framework. In particular, 
                we utilize the R/Shiny environment to minimize the loaded data. If you use 
                DNAlandscapeR in your research, please cite <REFERENCE?>...'),
        tags$h3('A WORD OF CAUTION-- The DNase data from ENCODE takes a while to load. 
                It will load though but slower than the other tracks. Be patient when
                selecting this. I will work to smooth the bigwig to make it faster.'),
        tags$h3('Descriptions of pre-loaded data:'),
        dataTableOutput('preloadedDataDescription'),
        tags$hr(), 
        width = 12
        )
)),
                   
tabPanel("Visualize", 
headerPanel(tags$h1(tags$b('Human Landscape'))),
  sidebarLayout(
  sidebarPanel(
    uiOutput("trackoptions"), 
    tags$hr(),
    textInput3("chr", HTML("<h5><b>Chr&nbsp;&nbsp;&nbsp;&nbsp;</b></h5>"), value = "12"),
    textInput3("start", HTML("<h5><b>Start&nbsp;</b></h5>"), value = "12660203"),
    textInput3("stop", HTML("<h5><b>Stop&nbsp;&nbsp;</b></h5>"), value = "12777226"),
    actionButton("plot.region", "Plot Region", style='padding:10px; font-size:80%'), 
    tags$hr(),
    textInput3("Gene", HTML("<h5><b>Gene&nbsp;</b></h5>"), value = "AGO3"),
    actionButton("plot.gene", "Plot Gene", style='padding:10px; font-size:80%'),
    tags$hr(),
    checkboxInput("showgenes", "Show Gene Annotation", value = TRUE, width = NULL),
    checkboxInput("showctcf", "Show CTCF Regions", value = FALSE, width = NULL),
    actionButton("zoom.in", "Zoom In", style='padding:10px; font-size:80%'),
    actionButton("zoom.out", "Zoom Out", style='padding:10px; font-size:80%'),
    actionButton("clear", "Clear", style='padding:10px; font-size:80%'),
    tags$hr(),
    actionButton("left.big", "<<", style='padding:10px; font-size:100%'),
    actionButton("left.small", "<", style='padding:10px; font-size:100%'),
    actionButton("right.small", ">", style='padding:10px; font-size:100%'),
    actionButton("right.big", ">>", style='padding:10px; font-size:100%'),
    tags$hr(),
    fileInput("skipRegions", "Upload Regions .bed File", multiple = FALSE, accept = NULL, width = NULL),
    conditionalPanel(condition = "input.skipRegions != NULL", textOutput("regionDescription")),
    actionButton("left.skip", "<<<", style='padding:10px; font-size:100%'),
    actionButton("right.skip", ">>>", style='padding:10px; font-size:100%'),
    tags$hr(),
    downloadButton("down", "Download Plot")
  ), 

  mainPanel(fluidRow(
    column(12, textOutput("regiontotal"), align = 'center')),
    plotOutput("plot", 
    dblclick = "plot_dblclick",
    brush = "plot_brush")),
  fluid = TRUE
)),

tabPanel("Upload",
pageWithSidebar(
    headerPanel(
        tags$h1(tags$b('Add local tracks'))
        ),
    sidebarPanel(
        fileInput("newTrack", h5(tags$b("Add new track:")), multiple = FALSE, accept = NULL, width = NULL),
        tags$hr(),
        textInput3("newTrackName", h5(tags$b("Specify track name:"), value = "")),
        tags$hr(),
        selectInput("datType", label = h5(tags$b("Data type:")), 
                choices = uploadchoices, #defined in global.R
                            selected = 1, width = "90%"),
        actionButton("addFile", "Add", style='padding:10px; font-size:80%')
        ),
    mainPanel(
        tags$h4('Successfuly loaded user files'),
        dataTableOutput('dt'),
        tags$hr()
        )
)),
theme = shinytheme("cosmo"),
footer = HTML('<P ALIGN=Center>DNAlandscapeR. &copy; <A HREF="mailto:caleblareau@g.harvard.edu">Caleb Lareau</A> & Martin Aryee'),
collapsible = TRUE, 
fluid = TRUE,
windowTitle = "DNAlandscapeR"
))
