source("helper.R")
source("global.R")

function(input, output, session) {
    val <- reactiveValues(region = NULL, curfil = NULL, alldat = NULL, list.tracks = f.list)
    observeEvent(input$plot.region, {
        val$region <- GRanges(seqnames=c(input$chr),ranges=IRanges(start=c(as.numeric(input$start)),end=c(as.numeric(input$stop))))
        })
    
    updateRegionVals <- function(){
        updateNumericInput(session, "chr", value = as.numeric(seqnames(val$region)))       
        updateNumericInput(session, "start", value = as.integer(start(ranges(range(val$region)))))  
        updateNumericInput(session, "stop", value = as.integer(end(ranges(range(val$region)))))     
    }
  
    observeEvent(input$plot.gene, {
        rda<-paste(system.file('rda',package='diffloop'),'human.genes.rda',sep='/')
        load(rda)
        t <- human.genes[mcols(human.genes)$id == as.character(input$Gene) ]
        val$region <- padGRanges(t, pad = as.integer(width(t)/2))
        updateRegionVals()
    })  
  
    observeEvent(input$zoom.out, {
        if (is.null(val$region)) return()
        val$region <- padGRanges(val$region, pad = as.integer(width(val$region)/2))
        updateRegionVals()
    })
    
     observeEvent(input$zoom.in, {
        if (is.null(val$region)) return()
        val$region <- padGRanges(val$region, pad = -1*as.integer(width(val$region)/4))
        updateRegionVals()
    })
     
      observeEvent(input$clear, {
        val$region <- NULL
    })
      
    observeEvent(input$down, {
        output$plotsave <- downloadHandler(
            filename = paste('plot-', Sys.Date(), '.pdf', sep=''),
            content = function(file){
                pdf(file = file, width=8.5, height=11)
                output$plot
                dev.off()
            })
    })
    
    p1 <- function(){  
        if (is.null(val$region)) return()
        if (length(input$tracks) == 0) return()
        par(mfrow=c(length(input$tracks)+input$showgenes, 1), oma = c(0, 0, 1, 0), mar = c(3, 5, 1, 1))
        for(i in input$tracks){
            i <- as.integer(i)
            if (i < 1000){ #ChIA-PET from Data Source File
                oneSampleLoopPlot(paste("data/loops/", c.full[[i]], sep = ""),val$region)
            } else if (i < 2000) { # BigWig Read Count Track
                bw.trackplot(paste("data/tracks/", e.full[[i-1000]], sep = ""), val$region)
            } else if (i < 3000){ #DNA Methylation
                methyl.bedgraph.trackplot(paste("data/methylation/", m.full[[i-2000]],sep = ""), val$region)
            } else {return()}
        }
        if(input$showgenes) humanAnnotation(val$region)
}
    output$down <- downloadHandler(
        filename <- function() {
            paste('plot-', Sys.Date(), '.pdf', sep='') },
        content <- function(file) {
            pdf(file, width = 8.5, height = 11)
            p1()
            dev.off()}
    )
      
    output$plot <- renderPlot({
        p1()
     }, height = 700)
    
    output$trackoptions <- renderUI({selectInput("tracks", label = h3("Tracks"), choices = val$list.tracks, selectize = TRUE, multiple= TRUE, selected = 0)})

    volumes <- getVolumes()
    shinyFileChoose(input, 'file', roots=volumes, session=session, restrictions=system.file(package='base'))
    output$filename <- renderPrint(as.character(parseFilePaths(volumes, input$file)$datapath))
    
    observeEvent(input$addFile, {
        #Error handling
        #validate(need(input$datType == "Loops" & input$fileformat == "rds" , "Loops object must be a .rds file"))
        
        #Update files
        val$curfil <- as.character(parseFilePaths(volumes, input$file)$datapath) 
        val$alldat <- as.matrix(rbind(val$alldat, cbind(val$curfil, input$fileformat, input$datType)), ncol = 3)
        colnames(val$alldat) <- c("File", "Format", "Data Type")
        dt <- as.data.frame(val$alldat)
        output$dt <- renderDataTable({dt})

        #Add track to global variables AND dynamic list option
        name <- basename(file_path_sans_ext(val$curfil))
        y <- val$list.tracks
        val <- 0
        if(input$datType == "Loops"){
            val <- as.list( max(unlist(y)[unlist(y) < 1000]) + 1 )
        } else if (input$datType == "Read.Depth"){
            val <- as.list( max(unlist(y)[unlist(y) < 2000]) + 1 )
        } else { #Methylation
            val <- as.list( max(unlist(y)[unlist(y) < 3000]) + 1 )
        }
        names(val) <- name
        val$list.tracks <- append(y, val)
        
        #Output new list
        output$trackoptions <- renderUI({selectInput("tracks", label = h3("Tracks"), choices = val$list.tracks, selectize = TRUE, multiple= TRUE, selected = 0)})

 })

    
}