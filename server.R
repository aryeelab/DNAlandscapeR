#Author: Caleb Lareau

source("plotter.R")
source("global.R")

function(input, output, session) {
    
    # Set up dataframe for data description tab
    dDF <- read.table("http://textuploader.com/5buvy/raw", header = TRUE, sep = "\t")
    dDF$PMID <- paste0('<a href="http://www.ncbi.nlm.nih.gov/pubmed/', dDF$PMID, '" target="_blank">', dDF$PMID, '</a>')
    output$preloadedDataDescription <- renderDataTable({dDF}, escape = FALSE)

    #---------------------------------#
    # Code for variable setup
    #---------------------------------#
    
    dynamic.val <- reactiveValues(
        region = GRanges(seqnames=c(default_chr),ranges=IRanges(start=c(default_start), end=c(default_end))),
        alldat = NULL,
        awsBuckets = NULL,

        # Initialize with human
        list.tracks = g_h.f.list, 
        c.full      = g_h.c.full,
        t.bw.full   = g_h.t.bw.full,
        t.bg.full   = g_h.t.bg.full,
        m.bw.full   = g_h.m.bw.full,
        m.bg.full   = g_h.m.bg.full, 
        i.full      = g_h.i.full,
        i.l.list    = NULL,
        c.list      = g_h.c.list,
        t.bw.list   = g_h.t.bw.list,
        t.bg.list   = g_h.t.bg.list,
        m.bw.list   = g_h.m.bw.list,
        m.bg.list   = g_h.m.bg.list,
        i.list      = g_h.i.list,
        i.res       = g_h.i.res,
        i.l.list    = NULL,
        organism    = "human",
        
        # Initialize human values
        h.f.list    = g_h.f.list, 
        h.c.full    = g_h.c.full,
        h.t.files   = g_h.t.files,
        h.m.files   = g_h.m.files,
        h.t.bw.full = g_h.t.bw.full,
        h.t.bg.full = g_h.t.bg.full,
        h.m.bw.full = g_h.m.bw.full,
        h.m.bg.full = g_h.m.bg.full,
        h.i.full    = g_h.i.full,
        h.i.l.full  = NULL,
        h.c.list    = g_h.c.list,
        h.t.bw.list = g_h.t.bw.list,
        h.t.bg.list = g_h.t.bg.list,
        h.m.bw.list = g_h.m.bw.list,
        h.m.bg.list = g_h.m.bg.list,
        h.i.list    = g_h.i.list,
        h.i.res     = g_h.i.res,
        h.i.l.list  = NULL,
        
        # Initialize mouse values
        m.f.list    = g_m.f.list,
        m.c.full    = g_m.c.full,
        m.t.files   = g_m.t.files,
        m.m.files   = g_m.m.files,
        m.t.bw.full = g_m.t.bw.full,
        m.t.bg.full = g_m.t.bg.full,
        m.m.bw.full = g_m.m.bw.full,
        m.m.bg.full = g_m.m.bg.full,
        m.i.full    = g_m.i.full,
        m.i.l.full  = NULL, 
        m.c.list    = g_m.c.list,
        m.t.bw.list = g_m.t.bw.list,
        m.t.bg.list = g_m.t.bg.list,
        m.m.bw.list = g_m.m.bw.list,
        m.m.bg.list = g_m.m.bg.list,
        m.i.list    = g_m.i.list,
        m.i.res     = g_m.i.res,
        m.i.l.list  = NULL,
        
        regionRow = 0,
        regions.df = NULL,
        fileAvail = FALSE,
        nregions = 0,
        track.names = NULL,
        start.tracks = 0,
        genes.avail = NULL
        )
    
    # Update current options to human
    updateToHuman <- function(){
        dynamic.val$organism <- "human"
        dynamic.val$list.tracks <- dynamic.val$h.f.list 
        dynamic.val$c.full <- dynamic.val$h.c.full
        dynamic.val$t.bw.full <- dynamic.val$h.t.bw.full
        dynamic.val$t.bg.full <- dynamic.val$h.t.bg.full
        dynamic.val$m.bw.full <- dynamic.val$h.m.bw.full
        dynamic.val$m.bg.full <- dynamic.val$h.m.bg.full
        dynamic.val$i.full <- dynamic.val$h.i.full
        dynamic.val$i.l.full <- dynamic.val$h.i.l.full
        dynamic.val$c.list <- dynamic.val$h.c.list
        dynamic.val$t.bw.list <- dynamic.val$h.t.bw.list
        dynamic.val$t.bg.list <- dynamic.val$h.t.bg.list
        dynamic.val$m.bw.list <- dynamic.val$h.m.bw.list
        dynamic.val$m.bg.list <- dynamic.val$h.m.bg.list
        dynamic.val$i.list <- dynamic.val$h.i.list  
        dynamic.val$i.l.list <- dynamic.val$h.i.l.list
        dynamic.val$i.res <- dynamic.val$h.i.res
    }
    
    updateToMouse <- function(){
        dynamic.val$organism <- "mouse"
        dynamic.val$list.tracks <- dynamic.val$m.f.list 
        dynamic.val$c.full <- dynamic.val$m.c.full
        dynamic.val$t.bw.full <- dynamic.val$m.t.bw.full
        dynamic.val$t.bg.full <- dynamic.val$m.t.bg.full
        dynamic.val$m.bw.full <- dynamic.val$m.m.bw.full
        dynamic.val$m.bg.full <- dynamic.val$m.m.bg.full
        dynamic.val$i.full    <- dynamic.val$m.i.full 
        dynamic.val$i.l.full <- dynamic.val$m.i.l.full
        dynamic.val$c.list <- dynamic.val$m.c.list
        dynamic.val$t.bw.list <- dynamic.val$m.t.bw.list
        dynamic.val$t.bg.list <- dynamic.val$m.t.bg.list
        dynamic.val$m.bw.list <- dynamic.val$m.m.bw.list
        dynamic.val$m.bg.list <- dynamic.val$m.m.bg.list
        dynamic.val$i.list <- dynamic.val$m.i.list 
        dynamic.val$i.l.list <- dynamic.val$m.i.l.list 
        dynamic.val$i.res <- dynamic.val$m.i.res
    }

    # Define dynamic variables based on organism selection
    observe({
        if(input$organism == 1)  updateToHuman()
        if(input$organism == 2)  updateToMouse()
    })
    
    output$trackoptions <- renderUI({selectInput("tracks", label = h3(tags$b("Select Tracks")),
                                                 choices = dynamic.val$list.tracks, selectize = TRUE,
                                                 multiple = TRUE, selected = dynamic.val$start.tracks)})
    
    output$plotGeneName <- renderUI({
        value <- "LMO2"
        if(input$organism == 2) value = "Grik4"
        textInput("Gene", label = HTML("<h4><b>Show Gene in Region </b></h4>"), value = value)
    })
    
    #Updates the text boxes of coordinates when a button is pressed.
    updateRegionVals <- function(){
        c <- as.character(data.frame(dynamic.val$region)[1,1])
        s <- as.numeric(start(ranges(range(dynamic.val$region))))
        e <- as.numeric(end(ranges(range(dynamic.val$region))))
        # Max chromosome sizes
        hg19.dist <- c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,
                  135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,
                  59128983,63025520,48129895,51304566,155270560)
        names(hg19.dist) <- c(paste0("chr", seq(1, 22, 1)), "chrX")
        
        mm9.dist <- c(197195432,181748087,159599783,155630120,152537259,149517037,152524553,131738871,124076172,
                      129993255,121843856,121257530,120284312,125194864,103494974,98319150,95272651,90772031,61342430,
                      166650296)
        names(mm9.dist) <- c(paste0("chr", seq(1, 19, 1)), "chrX")
        
        # Do some checking first to guard against 0, max, or switched input
        if(input$organism == 1){
            tchr <- c
            if(!(paste0("chr", c) %in% names(hg19.dist))){
                output$topText <- renderText(paste0("Chromosome ", tchr, " not found; plotting 1 instead"))
                c <- 1
            }
        } else {
            if(!(paste0("chr", c) %in% names(mm9.dist))){
                tchr <- c
                output$topText <- renderText(paste0("Chromosome ", tchr, " not found; plotting 1 instead"))
                c <- 1
            }
        }
        
        max_chrom_dist <- ifelse(input$organism == 1, hg19.dist[paste0("chr", c)], mm9.dist[paste0("chr", c)])
        s <- max(s, 0)
        e <- min(max_chrom_dist, e)
        
        dynamic.val$region <- GRanges(seqnames=c(c), ranges=IRanges(start=c(as.numeric(s)), end=c(as.numeric(e))))
        updateNumericInput(session, "chr", value = c)     
        updateNumericInput(session, "start", value = s)  
        updateNumericInput(session, "stop", value = e)
        updateTextInput(session, "ucscCoord", value = paste0("chr", c, ":", as.character(s), "-", as.character(e)))
    }
    
    observe({dynamic.val$track.names <- names(dynamic.val$list.tracks)[match(input$tracks, dynamic.val$list.tracks)]})

    
    #---------------------------------#
    # Code for various buttons
    #---------------------------------#
    
    observeEvent(input$initializeExample, {
        desired.tracks <- c("IMR90-HiC", "GM12878-HiC", "GM12878-ChIA-Pet-CTCF", "GM12878-ChIP-Seq-CTCF")
        dynamic.val$start.tracks <- g_h.f.list[desired.tracks]
    })
    
    observeEvent(input$plot.region, {
        dynamic.val$region <- GRanges(seqnames=c(input$chr),
                    ranges=IRanges(start=c(as.numeric(input$start)), end=c(as.numeric(input$stop))))
        makePlot()
        updateRegionVals()
    })
    
    observeEvent(input$plot.region2, {
        coord <- input$ucscCoord
        coord <- gsub("chr", "", coord)
        coord <- gsub(",", "", coord)
        coord <- gsub(":", " ", coord)
        coord <- gsub("-", " ", coord)
        coord <- gsub("\\t", " ", coord)
        coord <- gsub(" +", " ", coord)
        coord <- unlist(strsplit(coord, split=" "))
        c <- coord[1]
        s <- coord[2]
        e <- coord[3]
        if(s > e){ #switch the values
            ttmp <- e
            e <- s
            s <- ttmp
        }
        dynamic.val$region <- GRanges(seqnames=c(c),
                    ranges=IRanges(start=c(as.numeric(s)), end=c(as.numeric(e))))
        makePlot()
        updateRegionVals()
    })
    
    observeEvent(input$plot.gene, {
        if(input$organism == 1){ load("data/GenomeAnnotation/hg19/geneinfo.rda")}
        if(input$organism == 2){ load("data/GenomeAnnotation/mm9/geneinfo.rda")}
        t <- geneinfo[toupper(geneinfo$gene) == toupper(as.character(input$Gene)),]
        if(dim(t)[1] > 0){
            t.gr <- GRanges(t[c(1,2,3,4)])
            dynamic.val$region <- padGRanges(t.gr, pad = 150000)
            updateRegionVals()
            makePlot()
        } else {
            output$topText <- renderText("Gene name not found")
        }
    })  
  
    observeEvent(input$zoom.out, {
        if (is.null(dynamic.val$region)) return()
        dynamic.val$region <- padGRanges(dynamic.val$region, pad = as.integer(width(dynamic.val$region)/2))
        updateRegionVals()
        makePlot()
    })
    
     observeEvent(input$zoom.in, {
        if (is.null(dynamic.val$region)) return()
        dynamic.val$region <- padGRanges(dynamic.val$region, pad = -1*as.integer(width(dynamic.val$region)/4))
        updateRegionVals()
        makePlot()
    })
     
    observeEvent(input$left.big, {
        if (is.null(dynamic.val$region)) return()
        dynamic.val$region <- shift(dynamic.val$region, -9*width(dynamic.val$region)/10)
        updateRegionVals()
        makePlot()
    })

    observeEvent(input$left.small, {
        if (is.null(dynamic.val$region)) return()
        dynamic.val$region <- shift(dynamic.val$region, -3*width(dynamic.val$region)/10)
        updateRegionVals()
        makePlot()
    })

    observeEvent(input$right.small, {
        if (is.null(dynamic.val$region)) return()
        dynamic.val$region <- shift(dynamic.val$region, 3*width(dynamic.val$region)/10)
        updateRegionVals()
        makePlot()
    })

    observeEvent(input$right.big, {
        if (is.null(dynamic.val$region)) return()
        dynamic.val$region <- shift(dynamic.val$region, 9*width(dynamic.val$region)/10)
        updateRegionVals()
        makePlot()
    })
     
    observeEvent(input$clear, {
        dynamic.val$region <- NULL
        # add something here
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
            chr_temp <- as.character(as.matrix(data.frame(dynamic.val$region)[1,1]))
            dynamic.val$region <- GRanges(seqnames=c(chr_temp),
                                          ranges=IRanges(start=c(data.frame(dynamic.val$region)[1,2]+s),
                                          end=c(data.frame(dynamic.val$region)[1,2]+e)))
            updateRegionVals()
            makePlot()
        } else { return() }
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
        makePlot()
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
        makePlot()
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
        makePlot()
    })
    
    observeEvent( input$refresh, { makePlot() })
    observeEvent(input$refresh2, { makePlot() })
    

    p1 <- function(){  
        if (isolate(is.null(dynamic.val$region))) return()
        if (length(isolate(input$tracks)) == 0) return()
        sg <- ifelse(isolate(input$showgenes) == 0, 0, 1)
        par(mfrow=c(length(isolate(input$tracks)) + sg, 1),
                oma = c(0, 1, 3, 0), mar = c(3, 5, 3, 1))
        masterPlotter(isolate(input), isolate(dynamic.val))
    }
    
    ## Hi-C Resolution Display ##
    output$HiCresolutions <- renderUI({
        hictracks <- input$tracks[as.integer(input$tracks) <= 6000000 & as.integer(input$tracks) > 5000000] 
        plotSamples <- gsub("-HiC", "", names( dynamic.val$i.res[match(hictracks, dynamic.val$i.list)]))
        lapply(plotSamples, function(sample) {
                res <- sort(as.integer(dynamic.val$i.res[[sample]]))
                choices <- as.list(res)
                names(choices) <- res
                selectInput(paste0(sample, "HiCRes"), paste0('Specify ', sample, " Resolution"),
                    choices = choices, selected = min(res))
            })
    })
    
    ## Drop down buttons##
    output$flipMe <- renderUI({
        dropdownButton(label = "Flip Tracks", status = "default", width = 100,
        checkboxGroupInput(inputId = "flipper", label = "", choices = setNames(as.list(input$tracks), dynamic.val$track.names))
    )})
    
    output$showGenomeAnnotation <- renderUI({dropdownButton(label = "Genome Label", status = "default",
        checkboxGroupInput(inputId = "showGA", label = "", choices = setNames(as.list(input$tracks), dynamic.val$track.names),
        selected = setNames(as.list(input$tracks), dynamic.val$track.names) )
    )})
    
    ## Download Handlers ##
    
    output$downloadLoops <- downloadHandler(
        filename = function() { paste('DNAlandscapeR-loops-', Sys.Date(), '.tsv', sep='') },
        content = function(file) {
            write.table(getLoops(), file, quote = FALSE, sep = "\t", row.names = FALSE)
        }
    )
    
    output$downloadRDS <- downloadHandler(
        filename = function() { paste('DNAlandscapeR-allData-', Sys.Date(), '.rds', sep='') },
        content = function(file) {
            saveRDS(getAllData(), file)
        }
    )
    
    getLoops <- function(){
        masterPlotter(isolate(input), isolate(dynamic.val), loopsdl = TRUE)
    }
    
    getAllData <- function(){
        masterPlotter(isolate(input), isolate(dynamic.val), datadl = TRUE)
    }
    
    output$downQuick <- downloadHandler(
        filename <- function() {
            paste('DNAlandscapeR-plot-', Sys.Date(), '.pdf', sep='') },
        content <- function(file) {
            pdf(file, width = 8.5, height = 11)
            p1()
            dev.off()}
    )
    
    output$bm.down <- downloadHandler(
        filename <- function() {
            paste0("DNAlandscapeR-plot-", Sys.Date(), input$bm.type[1:3]) },
        content <- function(file) {
            bitmap(file, type = input$bm.type, height = as.numeric(input$bm.height), width = as.numeric(input$bm.width),
                   res = as.integer(input$bm.res), units = as.character(input$bm.units))
            p1()
            dev.off()}
    )
    
    makePlot <- function(){ 
        output$plot <- renderPlot({isolate(p1())}, height = 920)
        output$topText <- renderText("")
    }
    
    
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
            valu <- as.list(suppressWarnings(max(max(unlist(y)[unlist(y) < 2000000]), 1000000)) + 1 )
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
            valu <- as.list(suppressWarnings(max(max(unlist(y)[unlist(y) < 3000000]), 2000000)) + 1 )
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
            valu <- as.list(suppressWarnings(max(max(unlist(y)[unlist(y) < 4000000]), 3000000)) + 1 )
            names(valu) <- name
            if(input$organismUpload == 1){
                dynamic.val$h.m.bw.full <- c(dynamic.val$h.m.bw.full, curfile)
                dynamic.val$h.m.bw.list <- append(dynamic.val$h.m.bw.list, valu)
            } else {
                dynamic.val$m.m.bw.full <- c(dynamic.val$m.m.bw.full, curfile)
                dynamic.val$m.m.bw.list <- append(dynamic.val$m.m.bw.list, valu)
            }
            
        # Append bedgraph methylation file
        } else if (input$datType == 5) {
            valu <- as.list(suppressWarnings(max(max(unlist(y)[unlist(y) < 5000000]), 4000000)) + 1 )
            names(valu) <- name
            if(input$organismUpload == 1){
                dynamic.val$h.m.bg.full <- c(dynamic.val$h.m.bg.full, curfile)
                dynamic.val$h.m.bg.list <- append(dynamic.val$h.m.bg.list, valu)
            } else {
                dynamic.val$m.m.bg.full <- c(dynamic.val$m.m.bg.full, curfile)
                dynamic.val$m.m.bg.list <- append(dynamic.val$m.m.bg.list, valu)
            }
            
        # Hi C data; .l is for "local"
        } else { #datType is 6
            valu <- as.list(suppressWarnings(max(max(unlist(y)[unlist(y) < 7000000]), 6000000)) + 1 )
            names(valu) <- name
            if(input$organismUpload == 1){
                dynamic.val$h.i.l.full <- c(dynamic.val$h.i.l.full, curfile)
                dynamic.val$h.i.l.list <- append(dynamic.val$h.i.l.list, valu)
                dynamic.val$i.l.full <- dynamic.val$h.i.l.full
                dynamic.val$i.l.list <- dynamic.val$h.i.l.list
            } else {
                dynamic.val$m.i.l.full <- c(dynamic.val$m.i.l.full, curfile)
                dynamic.val$m.i.l.list <- append(dynamic.val$m.i.l.list, valu)
                dynamic.val$i.l.full <- dynamic.val$m.i.l.full
                dynamic.val$i.l.list <- dynamic.val$m.i.l.list
            }
        }
        
        # Update the tracks listed
        if(input$organismUpload == 1){
            dynamic.val$h.f.list <- append(dynamic.val$h.f.list, valu)
        } else {
            dynamic.val$m.f.list <- append(dynamic.val$m.f.list, valu)
        }        
    })
    
    observeEvent(input$addAWSBucket, {
        if(is.null(input$newBucket) | input$newBucket == "" ) return()
        bucket <-  gsub("^\\s+|\\s+$", "", input$newBucket) #gets rid of pesky white spaces
        ndf <- data.frame(paste0("http://s3.amazonaws.com/", bucket))
        names(ndf) <- "URL"
        dynamic.val$awsBuckets <- data.frame(rbind(dynamic.val$awsBuckets, ndf))
        names(dynamic.val$awsBuckets) <- "URL"
        dynamic.val <- importAmazonAWSBucket(bucket, dynamic.val)
        awsLoaded <- as.data.frame(dynamic.val$awsBuckets)
        output$awsLoaded <- renderDataTable({awsLoaded})
        
        # Display the correct values
        if(input$organism == 1){
            updateToHuman()
        } else {
            updateToMouse()
        }
        
        # Wipe record of individual local file upload
        dynamic.val$alldat <- NULL
        dt <- as.data.frame(dynamic.val$alldat)
        output$dt <- renderDataTable({dt})
        
    })
}
