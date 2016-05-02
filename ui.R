library(shiny)

pageWithSidebar(
  headerPanel("DNA Regional Viewer"),
  sidebarPanel(
    selectInput("tracks", label = h3("Tracks"), 
        choices = list(
            
            # Below 1,000-- ChIA Pet
            "jurkat_1-ChIA-PET" = 1,
            "jurkat_2-ChIA-PET" = 2,
            "naive_esc_1-ChIA-PET" = 3,
            "naive_esc_2-ChIA-PET" = 4,
            "primed_esc_1-ChIA-PET" = 5,
            "primed_esc_2-ChIA-PET" = 6,
            
            # From 1,000-2,000-- Epigenetics Tracks
            "Naive_1-H3k27ac" = 1001,
            "Primed_1-H3k27ac" = 1002,
            "Jurkat-H3k27ac" = 1003,
            
            # From 2,000-3,000-- DNA Methylation
            "Jurkat_450k" = 2001
            
            ), 
        selectize = TRUE, multiple= TRUE, selected = 0),
    textInput("chr", "Chromsome", value = "1", width = NULL, placeholder = NULL),
    textInput("start", "Start", value = "36100000", width = NULL, placeholder = NULL),
    textInput("stop", "Stop", value = "36300000", width = NULL, placeholder = NULL),
    textInput("Gene", "Gene", value = "AGO3", width = NULL, placeholder = NULL),
    actionButton("plot.region", "Plot Region"), 
    actionButton("plot.gene", "Plot Gene"),
    downloadButton("down", "Download Plot"),
    actionButton("clear", "Clear")
  ),

  mainPanel(
      plotOutput("plot")
      )
)

