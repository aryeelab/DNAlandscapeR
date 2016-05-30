#Author: Caleb Lareau

source("plotter.R")
source("global.R")

options(shiny.error=browser)
options(shiny.maxRequestSize=250*1024^2) #250 MB Max file size; can reparameterize if needed
options(warn=-1)

function(input, output, session) {
    
    # Set up dataframe
    dDF <- read.table("data/data-description.txt", header = TRUE, sep = "\t")
    output$preloadedDataDescription <- renderDataTable({dDF})
    output$regiontotal <- renderText("")
    
    #---------------------------------#
    # Code for variable setup
    #---------------------------------#
    
    dynamic.val <- reactiveValues(
        region = NULL,
        alldat = NULL,
        acceptedGenes = NULL, 
        
        # Initialize with human
        list.tracks = g_h.f.list, 
        c.full    = g_h.c.full,
        t.bw.full = g_h.t.bw.full,
        t.bg.full = g_h.t.bg.full,
        m.bw.full = g_h.m.bw.full,
        m.bg.full = g_h.m.bg.full, 
        c.list    = g_h.c.list,
        t.bw.list = g_h.t.bw.list,
        t.bg.list = g_h.t.bg.list,
        m.bw.list = g_h.m.bw.list,
        m.bg.list = g_h.m.bg.list,
        organism = "human",
        
        # Initialize human values
        h.f.list    = g_h.f.list, 
        h.c.full    = g_h.c.full,
        h.t.bw.full = g_h.t.bw.full,
        h.t.bg.full = g_h.t.bg.full,
        h.m.bw.full = g_h.m.bw.full,
        h.m.bg.full = g_h.m.bg.full, 
        h.c.list    = g_h.c.list,
        h.t.bw.list = g_h.t.bw.list,
        h.t.bg.list = g_h.t.bg.list,
        h.m.bw.list = g_h.m.bw.list,
        h.m.bg.list = g_h.m.bg.list,
        
        # Initialize mouse values
        h.f.list    = g_m.f.list,
        m.c.full    = g_m.c.full,
        m.t.bw.full = g_m.t.bw.full,
        m.t.bg.full = g_m.t.bg.full,
        m.m.bw.full = g_m.m.bw.full,
        m.m.bg.full = g_m.m.bg.full, 
        m.c.list    = g_m.c.list,
        m.t.bw.list = g_m.t.bw.list,
        m.t.bg.list = g_m.t.bg.list,
        m.m.bw.list = g_m.m.bw.list,
        m.m.bg.list = g_m.m.bg.list,
        
        regionRow = 0,
        regions.df = NULL,
        fileAvail = FALSE,
        nregions = 0
        )

    # Define dynamic variables based on organism selection
    observe({
    if(input$organism == 1) {
        dynamic.val$organism <- "human"
        dynamic.val$list.tracks <- dynamic.val$h.f.list 
        dynamic.val$c.full <- dynamic.val$h.c.full
        dynamic.val$t.bw.full <- dynamic.val$h.t.bw.full
        dynamic.val$t.bg.full <- dynamic.val$h.t.bg.full
        dynamic.val$m.bw.full <- dynamic.val$h.m.bw.full
        dynamic.val$m.bg.full <- dynamic.val$h.m.bg.full 
        dynamic.val$c.list <- dynamic.val$h.c.list
        dynamic.val$t.bw.list <- dynamic.val$h.t.bw.list
        dynamic.val$t.bg.list <- dynamic.val$h.t.bg.list
        dynamic.val$m.bw.list <- dynamic.val$h.m.bw.list
        dynamic.val$m.bg.list <- dynamic.val$h.m.bg.list
    }
    
    if(input$organism == 2) {
        dynamic.val$organism <- "mouse"
        dynamic.val$list.tracks <- dynamic.val$m.f.list 
        dynamic.val$c.full <- dynamic.val$m.c.full
        dynamic.val$t.bw.full <- dynamic.val$m.t.bw.full
        dynamic.val$t.bg.full <- dynamic.val$m.t.bg.full
        dynamic.val$m.bw.full <- dynamic.val$m.m.bw.full
        dynamic.val$m.bg.full <- dynamic.val$m.m.bg.full 
        dynamic.val$c.list <- dynamic.val$m.c.list
        dynamic.val$t.bw.list <- dynamic.val$m.t.bw.list
        dynamic.val$t.bg.list <- dynamic.val$m.t.bg.list
        dynamic.val$m.bw.list <- dynamic.val$m.m.bw.list
        dynamic.val$m.bg.list <- dynamic.val$m.m.bg.list
    }
    })
    
    #Updates the text boxes of coordinates when a button is pressed.
    updateRegionVals <- function(){
        updateNumericInput(session, "chr", value = data.frame(dynamic.val$region)[1,1])     
        updateNumericInput(session, "start", value = as.integer(start(ranges(range(dynamic.val$region)))))  
        updateNumericInput(session, "stop", value = as.integer(end(ranges(range(dynamic.val$region)))))   
        output$regiontotal <-renderText(paste(paste(paste("chr", input$chr, sep = ""),
                                        input$start, sep = ":"),input$stop, sep = "-"))
    }
    
    # Gets all gene names for the current region
    updateDisplayedGenes <- function(){
        if(input$showgenes){
            chrom <- as.character(seqnames(dynamic.val$region))
            chromchr <- paste(c("chr", as.character(chrom)), collapse = "")
            start <- as.integer(start(ranges(range(dynamic.val$region))))
            end <- as.integer(end(ranges(range(dynamic.val$region))))
            
            geneinfo <- data.frame()
            if(input$organism == 1) load("data/GenomeAnnotation/hg19/geneinfo.rda")
            if(input$organism == 2) load("data/GenomeAnnotation/mm9/geneinfo.rda")
            
            geneinfo <- geneinfo[geneinfo$chrom == chrom & geneinfo$start > start & geneinfo$stop < end,]
            dynamic.val$acceptedGenes <- unique(geneinfo$gene)
        }
    }
    
    #---------------------------------#
    # Code for various buttons
    #---------------------------------#
    
    observeEvent(input$plot.region, {
        dynamic.val$region <- GRanges(seqnames=c(input$chr),
                                      ranges=IRanges(start=c(as.numeric(input$start)),
                                                     end=c(as.numeric(input$stop))))
        output$regiontotal <-renderText(paste(paste(paste("chr", input$chr, sep = ""),
                                        input$start, sep = ":"),input$stop, sep = "-"))
    })
  
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
     
    observeEvent(input$left.big, {
        if (is.null(dynamic.val$region)) return()
        dynamic.val$region <- shift(dynamic.val$region, -9*width(dynamic.val$region)/10)
        updateRegionVals()
    })

    observeEvent(input$left.small, {
        if (is.null(dynamic.val$region)) return()
        dynamic.val$region <- shift(dynamic.val$region, -3*width(dynamic.val$region)/10)
        updateRegionVals()
    })

    observeEvent(input$right.small, {
        if (is.null(dynamic.val$region)) return()
        dynamic.val$region <- shift(dynamic.val$region, 3*width(dynamic.val$region)/10)
        updateRegionVals()
    })

    observeEvent(input$right.big, {
        if (is.null(dynamic.val$region)) return()
        dynamic.val$region <- shift(dynamic.val$region, 9*width(dynamic.val$region)/10)
        updateRegionVals()
    })
     
    observeEvent(input$clear, {
        dynamic.val$region <- NULL
    })
      
    observeEvent(input$plot_dblclick, {
        brush <- input$plot_brush
        if (!is.null(brush)) {
            fulldist <- brush$domain$right - brush$domain$left
            
            # Compute where the brush is occupying relative to window
            startprop <- (brush$xmin - brush$domain$left)/fulldist
            endprop <- (brush$xmax - brush$domain$left)/fulldist
            
            # Map that proportion to the GRanges width
            s <- as.integer(width(dynamic.val$region) * startprop)
            e <- as.integer(width(dynamic.val$region) * endprop)
            chr_temp <- as.numeric(as.matrix(data.frame(dynamic.val$region)[1,1]))
            dynamic.val$region <- GRanges(seqnames=c(chr_temp),
                                          ranges=IRanges(start=c(data.frame(dynamic.val$region)[1,2]+s),
                                          end=c(data.frame(dynamic.val$region)[1,2]+e)))
            updateRegionVals()
        } else { return() }
    })

    observeEvent(input$down, {
        output$plotsave <- downloadHandler(
            filename = paste('DNAlandscapeR-', Sys.Date(), '-plot.pdf', sep=''),
            content = function(file){
                pdf(file = file, width=8.5, height=11)
                output$plot
                dev.off()
            })
    })
    
    observeEvent(input$skipRegions, {
         if(is.null(input$skipRegions)){
                dynamic.val$fileAvail <- FALSE
            } else {
                dynamic.val$fileAvail <- TRUE
                dynamic.val$regions.df <- read.table(input$skipRegions$datapath)   
                dynamic.val$nregions <- dim(dynamic.val$regions.df)[1]
                output$regionDescription <- renderText(paste("Displaying region " , as.character(
                    dynamic.val$regionRow), " of ", as.character(dynamic.val$nregions), sep = ""))
            }
    })

    observeEvent(input$left.skip, {
        if (!dynamic.val$fileAvail | dynamic.val$regionRow <= 1) return()
        dynamic.val$regionRow <- dynamic.val$regionRow - 1
        chr <- dynamic.val$regions.df[dynamic.val$regionRow,1]
        chr <- gsub("^chr(.*)$", "\\1", chr)
        dynamic.val$region <- GRanges(seqnames=c(chr),
                                ranges=IRanges(start=c(dynamic.val$regions.df[dynamic.val$regionRow,2]),
                                end=c(dynamic.val$regions.df[dynamic.val$regionRow,3])))
        updateRegionVals()
    })
    
    observeEvent(input$right.skip, {
        if (!dynamic.val$fileAvail | dynamic.val$regionRow == dynamic.val$nregions) return()
        dynamic.val$regionRow <- dynamic.val$regionRow + 1
        chr <- dynamic.val$regions.df[dynamic.val$regionRow,1]
        chr <- gsub("^chr(.*)$", "\\1", chr)
        dynamic.val$region <- GRanges(seqnames=c(chr),
                                ranges=IRanges(start=c(dynamic.val$regions.df[dynamic.val$regionRow,2]),
                                end=c(dynamic.val$regions.df[dynamic.val$regionRow,3])))
        updateRegionVals()
    })
    
    observeEvent(input$showgenes, {
        if(!input$showgenes){
             dynamic.val$acceptedGenes <- NULL
        }
    })
    
    p1 <- function(){  
        if (is.null(dynamic.val$region)) return()
        if (length(input$tracks) == 0) return()
        updateDisplayedGenes()
        par(mfrow=c(length(input$tracks) + input$showgenes, 1),
                oma = c(0, 0, 1, 0), mar = c(3, 5, 1, 1))
        masterPlotter(input, dynamic.val)
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
    
    output$trackoptions <- renderUI({selectInput("tracks", label = h3(tags$b("Select Tracks")),
                                                 choices = dynamic.val$list.tracks, selectize = TRUE,
                                                 multiple = TRUE, selected = 0)})
    
    output$specifiedGenes <- renderUI({selectInput("plotGenes", label = h4(tags$b("Displayed Genes")),
                                                 choices = dynamic.val$acceptedGenes, selectize = TRUE,
                                                 multiple = TRUE, selected = dynamic.val$acceptedGenes)})

    #---------------------------------#
    # Code for input tab
    #---------------------------------#
    
    volumes <- getVolumes()
    shinyFileChoose(input, 'file', roots=volumes, session=session, restrictions=system.file(package='base'))
    output$filename <- renderPrint({
        fnc <- as.character(parseFilePaths(volumes, input$file)$datapath)
        if(identical(fnc, character(0))){ "" } else { fnc }
        })
    
    observeEvent(input$addFile, {
        if(is.null(input$newTrack)) return()
        curfile <- input$newTrack$datapath
        
        #Update data frame
        orgnsm <- "Human"
        if(input$organismUpload == 1){
            orgnsm <- "Human"
        } else {
            orgnsm <- "Mouse"  
        }
        dynamic.val$alldat <- as.matrix(rbind(dynamic.val$alldat,
                        cbind(input$newTrackName, orgnsm, names(uploadchoices)[as.numeric(input$datType)])), ncol = 2)
        colnames(dynamic.val$alldat) <- c("Sample", "Organism", "Data Type")
        dt <- as.data.frame(dynamic.val$alldat)
        output$dt <- renderDataTable({dt})

        #Add track to reactive variable for correct organism
        name <- input$newTrackName
        if(input$organismUpload == 1){
            y <- dynamic.val$h.f.list
        } else {
            y <- dynamic.val$m.f.list  
        }
        valu <- 0
        
        # Append ChIA-PET
        if(input$datType == 1){
            valu <- as.list(suppressWarnings(max(max(unlist(y)[unlist(y) < 1000000]), 0)) + 1 )
            names(valu) <- name
            if(input$organismUpload == 1){
                dynamic.val$h.c.full <- c(dynamic.val$h.c.full, curfile)
                dynamic.val$h.c.list <- append(dynamic.val$h.c.list, valu)
            } else {
                dynamic.val$m.c.full <- c(dynamic.val$m.c.full, curfile)
                dynamic.val$m.c.list <- append(dynamic.val$m.c.list, valu) 
            }
        
        # Append track/bigwig file        
        } else if (input$datType == 2) {
            valu <- as.list(suppressWarnings(max(max(unlist(y)[unlist(y) < 2000000]), 0)) + 1 )
            names(valu) <- name
            if(input$organismUpload == 1){
                dynamic.val$h.t.bw.full <- c(dynamic.val$h.t.bw.full, curfile)
                dynamic.val$h.t.bw.list <- append(dynamic.val$h.t.bw.list, valu)
            } else {
                dynamic.val$m.t.bw.full <- c(dynamic.val$m.t.bw.full, curfile)
                dynamic.val$m.t.bw.list <- append(dynamic.val$m.t.bw.list, valu)
            }
        
        # Append bedgraph track file
        } else if (input$datType == 3) {
            valu <- as.list(suppressWarnings(max(max(unlist(y)[unlist(y) < 3000000]), 0)) + 1 )
            names(valu) <- name
            if(input$organismUpload == 1){
                dynamic.val$h.t.bg.full <- c(dynamic.val$h.t.bg.full, curfile)
                dynamic.val$h.t.bg.list <- append(dynamic.val$h.t.bg.list, valu)
            } else {
                dynamic.val$m.t.bg.full <- c(dynamic.val$m.t.bg.full, curfile)
                dynamic.val$m.t.bg.list <- append(dynamic.val$m.t.bg.list, valu)
            }
        
        # Append bigwig methylation file
        } else if (input$datType == 4) {
            valu <- as.list(suppressWarnings(max(max(unlist(y)[unlist(y) < 4000000]), 0)) + 1 )
            names(valu) <- name
            if(input$organismUpload == 1){
                dynamic.val$h.m.bw.full <- c(dynamic.val$h.m.bw.full, curfile)
                dynamic.val$h.m.bw.list <- append(dynamic.val$h.m.bw.list, valu)
            } else {
                dynamic.val$m.m.bw.full <- c(dynamic.val$m.m.bw.full, curfile)
                dynamic.val$m.m.bw.list <- append(dynamic.val$m.m.bw.list, valu)
            }
            
        # Append bedgraph methylation file
        } else { 
            valu <- as.list(suppressWarnings(max(max(unlist(y)[unlist(y) < 5000000]), 0)) + 1 )
            names(valu) <- name
            if(input$organismUpload == 1){
                dynamic.val$h.m.bg.full <- c(dynamic.val$h.m.bg.full, curfile)
                dynamic.val$h.m.bg.list <- append(dynamic.val$h.m.bg.list, valu)
            } else {
                dynamic.val$m.m.bg.full <- c(dynamic.val$m.m.bg.full, curfile)
                dynamic.val$m.m.bg.list <- append(dynamic.val$m.m.bg.list, valu)
            }
        }
        
        # Update the tracks listed
        if(input$organismUpload == 1){
            dynamic.val$h.f.list <- append(dynamic.val$h.f.list, valu)
        } else {
            dynamic.val$m.f.list <- append(dynamic.val$m.f.list, valu)
        }        
    })
}
