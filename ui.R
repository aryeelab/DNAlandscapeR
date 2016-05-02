library(shiny)
chiapet_input <- list(
            
            # Below 1,000-- ChIA Pet
            "jurkat_1-ChIA-PET" = 1,
            "jurkat_2-ChIA-PET" = 2,
            "naive_esc_1-ChIA-PET" = 3,
            "naive_esc_2-ChIA-PET" = 4,
            "primed_esc_1-ChIA-PET" = 5,
            "primed_esc_2-ChIA-PET" = 6
            
            )
# From 1,000-2,000-- Epigenetics Tracks
e.full <- list.files("data/tracks/")
e.names <- file_path_sans_ext(e.full)
e.list <- as.list(seq(1,length(e.names), by = 1) + 1000)
names(e.list) <- e.names

# From 2,000-3,000-- DNA Methylation
m.full <- list.files("data/methylation/")
m.names <- file_path_sans_ext(m.full)
m.list <- as.list(seq(1,length(m.names), by= 1) + 2000)
names(m.list) <- m.names

f.list <- append(chiapet_input, append(e.list, m.list))

pageWithSidebar(
  headerPanel("DNA Regional Viewer"),
  sidebarPanel(
    selectInput("tracks", label = h3("Tracks"), 
        choices = f.list, 
        selectize = TRUE, multiple= TRUE, selected = 0),
    textInput("chr", "Chromsome", value = "1", width = NULL, placeholder = NULL),
    textInput("start", "Start", value = "36100000", width = NULL, placeholder = NULL),
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

