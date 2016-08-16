# DNAlandscapeR UI # 
source("global.R")

shinyUI(navbarPage(HTML("<img src='harvard-logo.png'/>"),
                   
tabPanel("Visualize", 
headerPanel(fluidRow(column(6, tags$h1(tags$b('DNA Landscape'))),
                     tags$h4(column(6, tags$h3(textOutput("topText")))))),
  sidebarLayout(
  sidebarPanel(
    uiOutput("trackoptions"), 
    tags$hr(),
    textInput("ucscCoord", HTML("<h4><b>Select Region Coordinates</b></h4>"), value = ucsc_coord),
    actionButton("plot.region2", "Plot Region", style='padding:10px; font-size:80%'),
    tags$br(),tags$hr(),
    uiOutput("plotGeneName"),
    actionButton("plot.gene", "Plot Gene", style='padding:10px; font-size:80%'),
    tags$br(),
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
    #tags$hr(),
    #fileInput("skipRegions", HTML("<h4><b>Jump through Regions</b></h4>"), multiple = FALSE, accept = NULL, width = NULL),
    #actionButton("left.skip", "<<<", style='padding:10px; font-size:100%'),
    #actionButton("right.skip", ">>>", style='padding:10px; font-size:100%'),
    tags$hr(),
    uiOutput("ucscGo"),
    tags$br(),
    downloadButton("downQuick", "Download Plot"), tags$br(), tags$br(),
    conditionalPanel(condition="input.initializeExample==0", uiOutput("initExamp")),
    tags$hr()), 

  mainPanel( plotOutput("plot", dblclick = "plot_dblclick",brush = "plot_brush") ),
  fluid = TRUE
),

bsCollapse(id = "collapseAdvancedPlotOptions", open = c("Panel1", "Panel2", "Panel3"), multiple = TRUE,
    bsCollapsePanel(title = HTML("<h4><b>Advanced Options</b></h4>"), value = "Panel1",
    fluidRow(
       column(4, radioButtons("organism", HTML("<h3><b>Specify Organism</b></h3>"),
                    choices = list("Human" = 1, "Mouse" = 2), selected = 1), tags$br(),
              selectInput("showgenes", label=HTML("<h4><b>Genome Annotation</b></h4>"),
                          choices = list("Gene Bodies" = 1, "Detailed Gene Annotation" = 2, "None" = 0), selected = 1)
              ),
       column(4,
              sliderInput("smoother", HTML("<h3><b>Smooth Epigenetic Peaks<br></b></h3>
                                           <h4><b>Window Size</b></h4>"), min=0, max=5000, value=500, step=50),
              selectInput("FUN", label=HTML("<h4><b>Function to Smooth</b></h4>"),
              choices = list("Mean"="mean", "Max"="max", "Median"="median"), selected = "mean"),
              checkboxInput("methylSmooth", "Also smooth WGBS Tracks", value = TRUE, width = NULL),
              checkboxInput("log2BW", "Log Transform Continuous Tracks", value = FALSE, width = NULL)
              ),
    column(4, HTML("<h3><b>Track Customization</b></h3>"),tags$br(),
              uiOutput("flipMe"), tags$br(),tags$br(),
              uiOutput("showGenomeAnnotation")
           )
     ),
    tags$hr(),
    fluidRow(
        column(4, textInput3("chr", HTML("<h5><b>Chr&nbsp;&nbsp;&nbsp;&nbsp;</b></h5>"), value = default_chr)),
        column(4, textInput3("start", HTML("<h5><b>Start&nbsp;</b></h5>"), value = default_start)),
        column(4, textInput3("stop", HTML("<h5><b>Stop&nbsp;&nbsp;</b></h5>"), value = default_end))),tags$hr(),
        actionButton("plot.region", "Plot Region"),
    style = "default"),
    
    bsCollapsePanel(title = HTML("<h4><b>ChIA-PET Tracks Configuration</b></h4>"), value = "Panel2",
    fluidRow(
        column(4,HTML("<h4><b>Individual Sample PET Cutoff</b></h4>"),
               uiOutput("petThresholds")),
        column(4, HTML("<h4><b>Configure ChIA-PET Visualization </b></h4>"),
              checkboxInput("colorLoops", "Color Loops Based on Biological Annotation", value = TRUE, width = NULL),
              checkboxInput("showSingleAnchors", "Show Single Anchors in ChIA-PET Tracks", value = FALSE, width = NULL)),
        column(4, selectInput("loopWidthNorm", label=HTML("<h4><b>ChIA-PET Loop Width Normalization</b></h4>"),
                           choices = list("Between Track" = 1, "Within Track" = 2, "None" = 0), selected = 1))),
    tags$hr(),
    actionButton("refresh2", "Refresh Plot"),
    style = "default"),
    
    bsCollapsePanel(title = HTML("<h4><b>Hi-C Tracks Configuration</b></h4>"), value = "Panel3",
    fluidRow(
        column(4,HTML("<h4><b>Individual Sample Resolution</b></h4>"),
               uiOutput("HiCresolutions")),
        column(4, HTML("<h4><b>Configure Hi-C Data</b></h4>"),
               checkboxInput("showlegend", "Show Legend on Plots", value = TRUE, width = NULL),
               radioButtons("HiCcutoff", HTML("<h4><b>Set Cutoff Thresholds</b></h4>"),
                    choices = list("Max/Min" = 1, "Quantile" = 2, "None" = 3), selected = 2),
    conditionalPanel(
            condition = "input.HiCcutoff == 1",
            textInput3("HiCmax", HTML("<h5><b>Max Value&nbsp;</b></h5>"), value = 100)
        ),
        conditionalPanel(
            condition = "input.HiCcutoff == 1",
            textInput3("HiCmin", HTML("<h5><b>Min Value&nbsp;&nbsp;</b></h5>"), value = 0)
        ),
        conditionalPanel(
            condition = "input.HiCcutoff == 2",
            textInput3("quantMax", HTML("<h5><b>Max Quant.&nbsp;</b></h5>"), value = 95)
        ),
        conditionalPanel(
            condition = "input.HiCcutoff == 2",
            textInput3("quantMin", HTML("<h5><b>Min Quant.&nbsp;&nbsp;</b></h5>"), value = 5)
        )),
        column(4, selectInput("HiCcolor",  HTML("<h4><b>Select Hi-C Color Theme</b></h4>"),
                    choices = color.choices, selected = 16),
                  selectInput("missingco",  HTML("<h4><b>Specify Missing Data Color</b></h4>"),
                    choices = missingco.choices, selected = "min"),
                  checkboxInput("log2hic", "Log Transform Hi-C Values", value = TRUE, width = NULL)
               )
        ),
    tags$hr(),
    actionButton("refresh3", "Refresh Plot"),
    style = "default"),
    
    bsCollapsePanel(title = HTML("<h4><b>Advanced Downloading</b></h4>"), value = "bmPanel",
    HTML('<h3><b><P ALIGN=Center>Plot Rendering</b></h3>'), tags$br(),
    fluidRow(
        column(1, tags$br()), 
        column(5,HTML('<h4><b>Parameters for the <a href="https://stat.ethz.ch/R-manual/R-devel/library/grDevices/html/dev2bitmap.html"
                      target="_blank">bitmap</a> function</b></h4>'),tags$br(),tags$br(),
                 textInput3("bm.type",HTML("<h5><b>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
                                           &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Type&nbsp;</b></h5>"), value = "png16m"),
                 textInput3("bm.res", HTML("<h5><b>Resolution&nbsp;</b></h5>"), value = 72), tags$br(),tags$br(),
                 downloadButton("bm.down", "Download Custom Plot")
               ),
        column(6,HTML('<h4><b>Plot Dimensions</b></h4>'),
                  textInput3("bm.height",HTML("<h5><b>Height&nbsp;</b></h5>"), value = 8.5),
                  textInput3("bm.width", HTML("<h5><b>&nbsp;&nbsp;Width&nbsp;</b></h5>"), value = 11),
                  textInput3("bm.units", HTML("<h5><b>&nbsp;&nbsp;&nbsp;Units&nbsp;</b></h5>"), value = "in")
               )),
    tags$hr(),tags$br(),
    HTML('<h3><b><P ALIGN=Center>Data Download</b></h3>'),
    fluidRow(
        column(1, tags$br()),
        column(10, HTML("<h5>The left button downloads a .tsv file of all loops that are currently displayed. 
                        The right button downloads a .rds of an R list that summarizes all data being displayed. </b>")),
        column(1, tags$br())),tags$br(),
    fluidRow(
        column(1, tags$br()),
        column(5, downloadButton("downloadLoops", "Download Displayed Loops")),
        column(5, downloadButton("downloadRDS", "Download All Displayed Data")),
        column(1, tags$br())),
    tags$hr(),
    style = "default")
)),

tabPanel("About",
fluidPage(
    headerPanel(
        HTML("<h2><b><P ALIGN=Center>Welcome to the DNAlandscapeR Epigenomics Browser</b></h2>
             <h3><b><P ALIGN=Center>Aryee Lab</b></h3>")
        ),
    mainPanel(includeHTML("www/welcome.html"),
              tags$br(),tags$br(),
              width = 12)
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

tabPanel("Import",
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
),
tags$hr(),
pageWithSidebar(
    headerPanel(
        tags$h1(tags$b('Import AWS Bucket'))
        ),
    sidebarPanel(
        textInput("newBucket", h5(tags$b("Specify bucket name:"), value = "")),
        tags$hr(),
        actionButton("addAWSBucket", "Import", style='padding:10px; font-size:80%')
        ),
    mainPanel(
        HTML("<h4><b>Loaded user buckets</b></h4>"),
        dataTableOutput('awsLoaded'),
        tags$hr()
        )
)


),
tabPanel("Guide",
    includeMarkdown("www/DNAlandscapeR-help.Rmd")
),
theme = shinytheme("cosmo"),
footer = HTML(paste0('<P ALIGN=Center>DNAlandscapeR. [Version <A HREF="https://github.com/aryeelab/DNAlandscapeR/tree/', sha, '">', short_sha, '</A>] &copy; <A HREF="mailto:caleblareau@g.harvard.edu">Caleb Lareau</A> & Martin Aryee')),
collapsible = TRUE, 
fluid = TRUE,
windowTitle = "DNAlandscapeR"
))
