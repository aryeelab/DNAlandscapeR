library(shiny)

pageWithSidebar(
  headerPanel("DNA Regional Viewer"),
  sidebarPanel(
    selectInput("tracks", label = h3("Tracks"), 
        choices = list("jurkat_1" = 1, "jurkat_2" = 2, "naive_esc_1" = 3, "naive_esc_2" = 4,
            "primed_esc_1" = 5, "primed_esc_2" = 6), 
        selectize = TRUE, multiple= TRUE, selected = 0),
    textInput("chr", "Chromsome", value = "1", width = NULL, placeholder = NULL),
    textInput("start", "Start", value = "36000000", width = NULL, placeholder = NULL),
    textInput("stop", "Stop", value = "36100000", width = NULL, placeholder = NULL),
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

