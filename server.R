source("plotter.R")
source("global.R")

function(input, output, session) {
    dynamic.val <- reactiveValues(
        region = NULL,
        curfil = NULL,
        alldat = NULL,
        list.tracks = f.list, 
        c.full = c.full,
        t.bw.full = t.bw.full,
        t.bg.full = t.bg.full,
        m.bw.full = m.bw.full,
        m.bg.full = m.bg.full
        )
    
    observeEvent(input$plot.region, {
        dynamic.val$region <- GRanges(seqnames=c(input$chr),ranges=IRanges(start=c(as.numeric(input$start)),end=c(as.numeric(input$stop))))
        })
    
    updateRegionVals <- function(){
        updateNumericInput(session, "chr", value = as.numeric(seqnames(dynamic.val$region)))       
        updateNumericInput(session, "start", value = as.integer(start(ranges(range(dynamic.val$region)))))  
        updateNumericInput(session, "stop", value = as.integer(end(ranges(range(dynamic.val$region)))))     
    }
  
    observeEvent(input$plot.gene, {
        rda<-paste(system.file('rda',package='diffloop'),'human.genes.rda',sep='/')
        load(rda)
        t <- human.genes[mcols(human.genes)$id == as.character(input$Gene) ]
        dynamic.val$region <- padGRanges(t, pad = as.integer(width(t)/2))
        updateRegionVals()
    })  
  
    observeEvent(input$zoom.out, {
        if (is.null(dynamic.val$region)) return()
        dynamic.val$region <- padGRanges(dynamic.val$region, pad = as.integer(width(dynamic.val$region)/2))
        updateRegionVals()
    })
    
     observeEvent(input$zoom.in, {
        if (is.null(dynamic.val$region)) return()
        dynamic.val$region <- padGRanges(dynamic.val$region, pad = -1*as.integer(width(dynamic.val$region)/4))
        updateRegionVals()
    })
     
      observeEvent(input$clear, {
        dynamic.val$region <- NULL
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
        if (is.null(dynamic.val$region)) return()
        if (length(input$tracks) == 0) return()
        par(mfrow=c(length(input$tracks)+input$showgenes, 1), oma = c(0, 0, 1, 0), mar = c(3, 5, 1, 1))
        for(i in input$tracks){
            i <- as.integer(i)
            if (i < 1000){ # ChIA-PET from rDS
                one.loopPlot(dynamic.val$c.full[[i]], dynamic.val$region)
            } else if (i < 2000) { # Track; BigWig
                bigwig.trackplot(dynamic.val$t.bw.full[[i-1000]], dynamic.val$region, "Read Depth")
            } else if (i < 3000){ # Track; Bedgraph
                bedgraph.trackplot(dynamic.val$t.bg.full[[i-2000]], dynamic.val$region, "Read Depth")
            } else if (i < 4000) { # Methyl; BigWig
                bigwig.trackplot(dynamic.val$m.bw.full[[i-3000]], dynamic.val$region, "Methylation")
            } else if (i < 5000){ # Methyl; Bedgraph
                bedgraph.trackplot(dynamic.val$m.bg.full[[i-4000]], dynamic.val$region, "Methylation")
            } else {return()}
        }
        if(input$showgenes) humanAnnotation(dynamic.val$region)
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
     }, height = 1000)
    
    output$trackoptions <- renderUI({selectInput("tracks", label = h3("Tracks"), choices = dynamic.val$list.tracks, selectize = TRUE, multiple = TRUE, selected = 0)})

    volumes <- getVolumes()
    shinyFileChoose(input, 'file', roots=volumes, session=session, restrictions=system.file(package='base'))
    output$filename <- renderPrint({
        fnc <- as.character(parseFilePaths(volumes, input$file)$datapath)
        if(identical(fnc, character(0))){ "" } else { fnc }
        })

    observeEvent(input$addFile, {
        #validate(need(if(input$datType == "Loops") input$fileformat == "rds" , "Loops input must be a .rds file"))
        #validate(need(if(input$datType == "Read.Depth") input$fileformat != "rds" , "Read.Depth cannot be a .rds file"))
        #validate(need(if(input$datType == "Methyl") input$fileformat != "rds" , "Methyl cannot be a .rds file"))

        #Update files
        dynamic.val$curfil <- as.character(parseFilePaths(volumes, input$file)$datapath) 
        dynamic.val$alldat <- as.matrix(rbind(dynamic.val$alldat, cbind(dynamic.val$curfil, input$fileformat, input$datType)), ncol = 3)
        colnames(dynamic.val$alldat) <- c("File", "Format", "Data Type")
        dt <- as.data.frame(dynamic.val$alldat)
        output$dt <- renderDataTable({dt})

        #Add track to global variables AND dynamic list option
        name <- basename(file_path_sans_ext(dynamic.val$curfil))
        y <- dynamic.val$list.tracks
        valu <- 0
        
        if(input$datType == "Loops"){
            valu <- as.list(suppressWarnings(max(max(unlist(y)[unlist(y) < 1000]), 0)) + 1 )
            names(valu) <- name
            dynamic.val$c.full <- c(dynamic.val$c.full, dynamic.val$curfil)
        } else if (input$datType == "Read.Depth" & input$fileformat == "BigWig") {
            valu <- as.list(suppressWarnings(max(max(unlist(y)[unlist(y) < 2000]), 0)) + 1 )
            names(valu) <- name
            dynamic.val$t.bw.full <- c(dynamic.val$t.bw.full, dynamic.val$curfil)
        } else if (input$datType == "Read.Depth" & input$fileformat == "Bedgraph") {
            valu <- as.list(suppressWarnings(max(max(unlist(y)[unlist(y) < 3000]), 0)) + 1 )
            names(valu) <- name
            dynamic.val$t.bg.full <- c(dynamic.val$t.bg.full, dynamic.val$curfil)
        } else if (input$datType == "Methyl" & input$fileformat == "BigWig") {
            valu <- as.list(suppressWarnings(max(max(unlist(y)[unlist(y) < 4000]), 0)) + 1 )
            names(valu) <- name
            dynamic.val$m.bw.full <- c(dynamic.val$m.bw.full, dynamic.val$curfil)
        } else { #Methyl; Bedgraph
            valu <- as.list(suppressWarnings(max(max(unlist(y)[unlist(y) < 5000]), 0)) + 1 )
            names(valu) <- name
            dynamic.val$m.full <- c(dynamic.val$m.full, dynamic.val$curfil)
        }

        #Output new list
        dynamic.val$list.tracks <- append(y, valu)
        output$trackoptions <- renderUI({selectInput("tracks", label = h3("Tracks"), choices = dynamic.val$list.tracks, selectize = TRUE, multiple = TRUE, selected = 0)})

    })

    
}