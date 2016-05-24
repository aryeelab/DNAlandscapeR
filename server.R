source("plotter.R")
source("global.R")

options(shiny.error=browser)
options(shiny.maxRequestSize=250*1024^2) #250 MB Max filze size
options(warn=-1)

function(input, output, session) {
    
    # Initialize Data Frame
    dDF <- read.table("data/data-description.txt", header = TRUE, sep = "\t")
    output$preloadedDataDescription <- renderDataTable({dDF})
    output$regiontotal <- renderText("")
    
    #Initialize Reactive variables
    dynamic.val <- reactiveValues(
        region = NULL,
        alldat = NULL,
        list.tracks = f.list, 
        c.full    = c.full,
        t.bw.full = t.bw.full,
        t.bg.full = t.bg.full,
        m.bw.full = m.bw.full,
        m.bg.full = m.bg.full, 
        c.list    = c.list,
        t.bw.list = t.bw.list,
        t.bg.list = t.bg.list,
        m.bw.list = m.bw.list,
        m.bg.list = m.bg.list, 
        regionRow = 0,
        regions.df = NULL,
        fileAvail = FALSE,
        nregions = 0
        )

    observeEvent(input$plot.region, {
        dynamic.val$region <- GRanges(seqnames=c(input$chr),
                                      ranges=IRanges(start=c(as.numeric(input$start)),
                                                     end=c(as.numeric(input$stop))))
        output$regiontotal <-renderText(paste(paste(paste("chr", input$chr, sep = ""),
                                        input$start, sep = ":"),input$stop, sep = "-"))
        })
    
    updateRegionVals <- function(){
        updateNumericInput(session, "chr", value = data.frame(dynamic.val$region)[1,1])     
        updateNumericInput(session, "start", value = as.integer(start(ranges(range(dynamic.val$region)))))  
        updateNumericInput(session, "stop", value = as.integer(end(ranges(range(dynamic.val$region)))))   
        output$regiontotal <-renderText(paste(paste(paste("chr", input$chr, sep = ""),
                                        input$start, sep = ":"),input$stop, sep = "-"))
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
        
    p1 <- function(){  
        if (is.null(dynamic.val$region)) return()
        if (length(input$tracks) == 0) return()
        par(mfrow=c(length(input$tracks) + input$showgenes + input$showctcf, 1),
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
     }, height = 800)
    
    output$trackoptions <- renderUI({selectInput("tracks", label = h3(tags$b("Tracks")),
                                                 choices = dynamic.val$list.tracks, selectize = TRUE,
                                                 multiple = TRUE, selected = 0)})

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
        dynamic.val$alldat <- as.matrix(rbind(dynamic.val$alldat,
                        cbind(input$newTrackName,names(uploadchoices)[as.numeric(input$datType)])), ncol = 2)
        colnames(dynamic.val$alldat) <- c("Sample", "Data Type")
        dt <- as.data.frame(dynamic.val$alldat)
        output$dt <- renderDataTable({dt})

        #Add track to global variables AND dynamic list option
        name <- input$newTrackName
        y <- dynamic.val$list.tracks
        valu <- 0
        
        if(input$datType == 1){
            valu <- as.list(suppressWarnings(max(max(unlist(y)[unlist(y) < 1000000]), 0)) + 1 )
            names(valu) <- name
            dynamic.val$c.full <- c(dynamic.val$c.full, curfile)
            dynamic.val$c.list <- append(dynamic.val$c.list, valu)
        } else if (input$datType == 2) {
            valu <- as.list(suppressWarnings(max(max(unlist(y)[unlist(y) < 2000000]), 0)) + 1 )
            names(valu) <- name
            dynamic.val$t.bw.full <- c(dynamic.val$t.bw.full, curfile)
            dynamic.val$t.bw.list <- append(dynamic.val$t.bw.list, valu)
        } else if (input$datType == 3) {
            valu <- as.list(suppressWarnings(max(max(unlist(y)[unlist(y) < 3000000]), 0)) + 1 )
            names(valu) <- name
            dynamic.val$t.bg.full <- c(dynamic.val$t.bg.full, curfile)
            dynamic.val$t.bg.list <- append(dynamic.val$t.bg.list, valu)
        } else if (input$datType == 4) {
            valu <- as.list(suppressWarnings(max(max(unlist(y)[unlist(y) < 4000000]), 0)) + 1 )
            names(valu) <- name
            dynamic.val$m.bw.full <- c(dynamic.val$m.bw.full, curfile)
            dynamic.val$m.bw.list <- append(dynamic.val$m.bw.list, valu)
        } else { #Methyl; Bedgraph; == 5
            valu <- as.list(suppressWarnings(max(max(unlist(y)[unlist(y) < 5000000]), 0)) + 1 )
            names(valu) <- name
            dynamic.val$m.bg.full <- c(dynamic.val$m.bg.full, curfile)
            dynamic.val$m.bg.list <- append(dynamic.val$m.bg.list, valu)
        }

        #Output new list
        dynamic.val$list.tracks <- append(y, valu)
        output$trackoptions <- renderUI({selectInput("tracks", label = h3(tags$b("Tracks")),
                choices = dynamic.val$list.tracks, selectize = TRUE, multiple = TRUE, selected = 0)})
    })
}