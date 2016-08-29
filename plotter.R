# Main functions for plotting the various tracks.
#region <- GRanges(seqnames=c("9"),ranges=IRanges(start=c(20912689),end=c(22216233)))


# masterPlotter is called from the server.R script and performs these operations:
# 1) subsets all ChIA-PET objects to determine the max counts for normalization
# 2) plots tracks calling later functions based on the indices.

# Additional parameters for download handling hijack this function and then exit
# before plotting to create data objects for downloading

readCachedRDS <- function(file) {
  if(grepl("amazonaws", file)){
    cached_file <- gsub("/", "__", sub(".*//", "", file))
    cached_file <- file.path(cache_dir, cached_file)
    if (file.exists(cached_file)) {
      x <- readRDS(cached_file)
      print(paste("Reading file from cache: ", cached_file))
    } else {
      filegz <- gzcon(url(file))
      x <- readRDS(filegz)
      close(filegz)
      saveRDS(x, file=cached_file)
      print(paste("Writing file to cache: ", cached_file))
    }
  } else { 
    x <- readRDS(file)
  }
  return(x)
}


masterPlotter <- function(input, dynamic.val, loopsdl = FALSE, datadl = FALSE){
    chia_pet_samples <- list() #Tracks samples linking with i
    chia_pet_objects <- c() #Tracks subsetted objects
    j <- 1 #Index of subsetted objects vector
    map_chia_pet.indices <- c() #maps i to order in vector of subsetted objects
    mc <- 1 #max counts
    one_anchor_samples <- list() #Tracks samples linking with i
    
    #Handle download if requested
    if(loopsdl) loopsTotal <- data.frame()
    if(datadl) datOut <- setNames(list(data.frame(dynamic.val$region)), "region")
    
    #First loop initalizes the ChIA-PET max values
    for(i in input$tracks){
        i <- as.integer(i)
        if (i < 1000000){ # ChIA-PET from RDS
            
            #Import object and subset
            file.conn <- dynamic.val$c.full[[i]]
            x <- readCachedRDS(file.conn)
            sample <- names(dynamic.val$c.list)[i]
            objReg <- removeSelfLoops(.subsetRegion.quick(x, dynamic.val$region, nanchors = 2))
            
            #Handle loops download if requested
            if(loopsdl){
                sdf <- summary(objReg)
                sdf$sample <- sample
                colnames(sdf)[7] <- "counts"
                loopsTotal <- rbind(loopsTotal, sdf)
            }
            
            #Grab loops with one anchor, in needed
            if(input$showSingleAnchors){
                oneAnchor <- .subsetRegion.quick(x, dynamic.val$region, nanchors = 1)
                oneAnchors <- oneAnchor@anchors[findOverlaps(oneAnchor@anchors, dynamic.val$region)@from]
                one_anchor_samples <- c(one_anchor_samples, oneAnchors)
            }
            
            #Update Max Counts
            if(max(objReg@counts) > mc) mc <- max(objReg@counts)
            
            #Add sample name to list
            valu_sample <- sample
            names(valu_sample) <- as.character(i)
            chia_pet_samples <- c(chia_pet_samples, valu_sample)
            
            #Add subsetted object to list
            chia_pet_objects <- append(chia_pet_objects, objReg)
            map_chia_pet.indices[j] <- i
            j <- j + 1
        }
    }
    
    # Change max counts if user wants
    if(input$loopWidthNorm == 0) mc = -3
    if(input$loopWidthNorm == 2) mc = -2
    
    if(loopsdl) return(loopsTotal)
    
    #Second loop does all the plotting
    for(i in input$tracks){ 
        flipped <- i %in% input$flipper
        showGA <- i %in% input$showGA
        i <- as.integer(i)
        if (i < 1000000) {
            sample <- chia_pet_samples[[as.character(i)]]
            minPets <- as.integer(input[[paste0(sample, "petThresh")]])
            oa <- try(one_anchor_samples[[which(map_chia_pet.indices == i)]], silent = TRUE)
            if("try-error" %in% class(oa)) oa <- NULL
            o <- one.loopPlot(objReg = chia_pet_objects[[which(map_chia_pet.indices == i)]], y = dynamic.val$region,
                         sample = sample, max_counts = mc, colorLoops = input$colorLoops,
                         oneAnchor = oa, flip = flipped, minPets = minPets, 
                         showGA = showGA, datadl = datadl)
            if(datadl) datOut <- append(datOut, setNames(list(o), chia_pet_samples[[as.character(i)]]))
        } else if (i < 2000000) { # Track; BigWig
            t <- i - 1000000
            sample <- names(dynamic.val$t.bw.list)[t]
            o <- bigwig.trackplot(dynamic.val$t.bw.full[[t]], dynamic.val$region, input$smoother, datadl = datadl,
                             FUN = input$FUN, "Depth", sample = sample, log2 = input$log2BW, flip = flipped, showGA = showGA)
            if(datadl) datOut <- append(datOut, setNames(list(o), sample))
        } else if (i < 3000000){ # Track; Bedgraph
            t <- i - 2000000
            sample <- names(dynamic.val$t.bg.list)[t]
            o <- bedgraph.trackplot(dynamic.val$t.bg.full[[t]], dynamic.val$region, "Depth", sample = sample, flip = flipped,
                               showGA = showGA, datadl = datadl)
            if(datadl) datOut <- append(datOut, setNames(list(o), sample))
        } else if (i < 4000000) { # Methyl; BigWig
            t <- i - 3000000
            sample <- names(dynamic.val$m.bw.list)[t]
            o <- bigwig.bumpPlot(dynamic.val$m.bw.full[[t]], dynamic.val$region, sample = sample, showGA = showGA, 
                            smoother = input$smoother, FUN = input$FUN, smoothBool = input$methylSmooth, flip = flipped,
                            datadl = datadl)
            if(datadl) datOut <- append(datOut, setNames(list(o), sample))
        } else if (i < 5000000){ # Methyl; Bedgraph
            t <- i - 4000000
            sample <- names(dynamic.val$m.bg.list)[t]
            o <- bedgraph.trackplot(dynamic.val$m.bg.full[[t]], dynamic.val$region, "Methylation", sample = sample, flip = flipped,
                     showGA = showGA,  datadl = datadl)
            if(datadl) datOut <- append(datOut, setNames(list(o), sample))
        } else if (i < 6000000){ # Hi-C Plot for Stuff that's preloaded on the server
            t <- i - 5000000
            sample.hic <- names(dynamic.val$i.list)[t]
            sample <- sample.hic
            fs <- dynamic.val$i.full
            res <- as.character(input[[paste0(sample, "HiCRes")]])
            chrom <- paste0("chr", as.character(seqnames(dynamic.val$region)))
            file <- fs[grepl(".rds", fs) & grepl(sample, fs)]
            hics4 <- readCachedRDS(file)
            hicdata <- hics4@resolutionNamedList[[res]][[chrom]]
            o <- hic.plot(hicdata, dynamic.val$region, sample = sample.hic, color = input$HiCcolor, log2trans = input$log2hic, flip = flipped,
                     missingco = input$missingco, showlegend = input$showlegend, showGA = showGA,  datadl = datadl, HiCmin = input$HiCmin,
                     HiCmax = input$HiCmax, custMaxMin = input$HiCcutoff, Qmin = input$quantMin, Qmax = input$quantMax)
            if(datadl) datOut <- append(datOut, setNames(list(o), sample.hic))
        } else if (i < 7000000) { # Local Hi-C plot    
            t <- i - 6000000
            chrom <- paste0("chr", as.character(seqnames(dynamic.val$region)))
            list.dat <- readRDS(dynamic.val$i.l.full[t])
            hicdata <- list.dat[[chrom]]
            sample <- names(dynamic.val$i.l.list)[t]
            o <- hic.plot(hicdata, dynamic.val$region, sample = sample, color = input$HiCcolor, log2trans = input$log2hic, flip = flipped,
                     missingco = input$missingco, showlegend = input$showlegend, showGA = showGA, datadl = datadl, HiCmin = input$HiCmin,
                     HiCmax = input$HiCmax, custMaxMin = input$HiCcutoff, Qmin = input$quantMin, Qmax = input$quantMax)
            if(datadl) datOut <- append(datOut, setNames(list(o), sample.hic))
        } else {return()}
    }
    e <- ifelse(input$showgenes == 2, TRUE, FALSE)
    if(input$showgenes > 0 & input$organism == 1){
        o <- geneAnnotation(dynamic.val$region, "human", exons = e, datadl)
        if(datadl) datOut <- append(datOut, setNames(list(o), "annotation"))
    }
    if(input$showgenes > 0 & input$organism == 2){
       o <-  geneAnnotation(dynamic.val$region, "mouse", exons = e, datadl)
       if(datadl) datOut <- append(datOut, setNames(list(o), "annotation"))
    }
    if(datadl) return(datOut)
}

# one.loopPlot has some specialized features for plotting only one sample's loops in these plots. The object is a loops object
# without self loops generated from the master function to determine the max_counts
one.loopPlot <- function(objReg, y, sample, max_counts, colorLoops = TRUE, oneAnchor = NULL,
                         flip, minPets, showGA, datadl) {
    
    if(datadl) return(summary(objReg))
    
    # Grab Regional Coordinates
    chrom <- as.character(seqnames(y))
    chromchr <- paste(c("chr", as.character(chrom)), collapse = "")
    start <- as.integer(start(ranges(range(y))))
    end <- as.integer(end(ranges(range(y))))
    
    # Make sure loop object is non-empty
    if(dim(objReg)[2] != 0){
        res <- objReg@rowData
        n <- dim(objReg@interactions)[1]  #number of interactions
        
        # Setup colors for plotting
        cs <- 0
        if(!is.null(res$loop.type) & colorLoops){
            cs <- res$loop.type
            cs <- gsub("e-p", "red", cs)
            cs <- gsub("p-p", "orange", cs)
            cs <- gsub("e-e", "mediumpurple1", cs)
            cs <- gsub("ctcf", "blue", cs)
            cs <- gsub("none", "black", cs)
        } else {
            cs <- rep("black", n)
        }
        
        # Setup Dataframe for Plot
        leftAnchor <- as.data.frame(objReg@anchors[objReg@interactions[,1]])[c(1, 2, 3)]
        LA <- do.call("rbind", replicate(1, leftAnchor, simplify = FALSE))
        rightAnchor <- as.data.frame(objReg@anchors[objReg@interactions[,2]])[c(1, 2, 3)]
        RA <- do.call("rbind", replicate(1, rightAnchor, simplify = FALSE))
        colnames(LA) <- c("chr_1", "start_1", "end_1")
        colnames(RA) <- c("chr_2", "start_2", "end_2")
        name <- rep(NA, n)
        strand_1 <- rep(".", n * 1)
        strand_2 <- rep(".", n * 1)
        score <- matrix(objReg@counts, ncol = 1)
        score[score[,1] < minPets ] <- 0
        bedPE <- data.frame(LA, RA, name, score, strand_1, strand_2, sample)

        w <- loopWidth(objReg)
        h <- sqrt(w/max(w))
        lwd <- 5 * (bedPE$score/max_counts)
        if(max_counts == -2){ #within track normalization
            max_counts <- max(bedPE$score)
            lwd <- 5 * (bedPE$score/max_counts)
        } else if (max_counts == -3) {
            max_counts <- 1
            lwd <- 3
        }

        # Add single loops
        if(!is.null(oneAnchor) ){
          if(dim(data.frame(oneAnchor))[1] != 0){
            #Make new data frame
            tdf <- data.frame(oneAnchor)
            a1df <- data.frame(
                chr_1 = tdf$seqnames, 
                start_1 = tdf$start,
                end_1 = tdf$start,
                chr_2 = tdf$seqnames,
                start_2 = tdf$end,
                end_2 = tdf$end,
                name = NA,
                score = max_counts,
                strand_1 = ".",
                strand_2 = ".",
                sample = sample
            )
            bedPE <- rbind(bedPE, a1df)
            
            #Update vectors
            cs <- c(cs, rep("forestgreen", dim(a1df)[1]))
            h <- c(h, rep(0.01, dim(a1df)[1]))
            lwd <- c(lwd, rep(4, dim(a1df)[1]))
          }
        }
        
        loplot <- recordPlot()
        print(bedPE)
        print(cs)
        plotBedpe(bedPE, chrom, start, end, color = cs, lwd = lwd, 
                  plottype = "loops", heights = h, lwdrange = c(0, 5), 
                  main = sample, adj=0, flip = flip)
        if(showGA) labelgenome(chromchr, start, end, side = 1, scipen = 20, n = 3, scale = "Mb", line = 0.18, chromline = 0.5, scaleline = 0.5)
        return(loplot)
        
    } else {
        # Return dummy plot
        loplot <- recordPlot()
        plotBedpe(data.frame(), chrom, start, end, color = c("blue"), lwd = 0, 
                  plottype = "loops", heights = 0, lwdrange = c(0, 0), 
                  main = sample, adj=0)
        if(showGA) labelgenome(chromchr, start, end, side = 1, scipen = 20, n = 3, scale = "Mb", line = 0.18, chromline = 0.5, scaleline = 0.5)
        return(loplot)   
    }
    
}

# bigwig.bumpPlot is used for methylation
bigwig.bumpPlot <- function(file, region, shade = TRUE, sample, showGA, smoother, FUN, smoothBool, flip, datadl){
    region.bed <- import.bw(file, which = addchr(region))
    
    if(smoothBool){
        tile <- unlist(tile(addchr(region), width = smoother)) 
        ovl <- findOverlaps(tile, region.bed)
        qh <- queryHits(ovl) 
        sh <- subjectHits(ovl) 
        values.t <- as.data.frame(tapply(mcols(region.bed[sh])$score, qh, get(FUN)))
        
        #A lot of extra effort to handle regions with no values
        colnames(values.t) <- "bwvalues"
        vNA <- data.frame(matrix(NA, ncol = 1, nrow = length(ranges(tile))))
        colnames(vNA) <- "NAss"
        ugly <- merge(vNA, values.t, by=0, all = TRUE, sort = F)
        ugly <- ugly[order(as.numeric(ugly$Row.names)), ]
        mcols(tile)$score <- unname(ugly$bwvalues, force = TRUE)
        region.bed <- suppressWarnings(tile[!is.na(mcols(tile)$score)])
    }
    
    region.bedgraph <- data.frame(region.bed)
    region.bedgraph <- region.bedgraph[,c(-4,-5)]
    if(datadl) return(region.bedgraph)
    if(flip) region.bedgraph$score <- region.bedgraph$score*(-1)

    chrom <- as.character(seqnames(region))
    chromchr <- paste(c("chr", as.character(chrom)), collapse = "")
    start <- as.integer(start(ranges(range(region))))
    end <- as.integer(end(ranges(range(region))))

    bumpplot <- recordPlot()

    pos <- region.bedgraph$start
    y <- region.bedgraph[,4]
    
    
    if(dim(region.bedgraph)[1] == 0) { # dummyplot
        plotBedpe(data.frame(), chrom, start, end, color = c("blue"), lwd = 0, plottype = "loops", heights = 0, lwdrange = c(0, 0), main = sample, adj=0)
    } else { #real plot
        if(!smoothBool){
            cluster_id <- clusterMaker(chr=chrom, pos=pos, maxGap = 100)
            smooth <- locfitByCluster(x=pos, y=y, cluster=cluster_id, bpSpan = 50)
            plot(pos, smooth$fitted, type="l", xaxt='n',bty = "n",xaxs="i",yaxs="i",main=sample,adj=0,ylab="")
        } else{
            plot(pos,y, type="l", xaxt='n',bty = "n",xaxs="i",yaxs="i",main=sample,adj=0,ylab="")
        }
    }
    
    if(showGA) labelgenome(chromchr, start, end, side = 1, scipen = 20, n = 3, scale = "Mb", line = 0.18, chromline = 0.5, scaleline = 0.5)
    if(shade & !flip) polygon(cbind(c(min(pos), pos, max(pos)), c(min(y), y, min(y))), border=NA, col="black")
    if(shade & flip) polygon(cbind(c(min(pos), pos, max(pos)), c(max(y), y, max(y))), border=NA, col="black")
    return(bumpplot)
}

# bigwig.trackplot is used for most epigenetic peaks
bigwig.trackplot <- function(file, region, smoother, datadl, FUN, ylab, sample, log2, flip, showGA){
    region.bed <- import.bw(file, which = addchr(region))
    # smooth
    if(smoother != 0){
        tile <- unlist(tile(addchr(region), width = smoother)) 
        ovl <- findOverlaps(tile, region.bed)
        qh <- queryHits(ovl) 
        sh <- subjectHits(ovl) 
        values.t <- as.data.frame(tapply(mcols(region.bed[sh])$score, qh, get(FUN)))
        
        #A lot of extra effort to handle regions with no values
        colnames(values.t) <- "bwvalues"
        vNA <- data.frame(matrix(NA, ncol = 1, nrow = length(ranges(tile))))
        colnames(vNA) <- "NAss"
        ugly <- merge(vNA, values.t, by=0, all = TRUE, sort = F)
        ugly <- ugly[order(as.numeric(ugly$Row.names)), ]
        mcols(tile)$score <- unname(ugly$bwvalues, force = TRUE)
        region.bed <- suppressWarnings(tile[!is.na(mcols(tile)$score)])
    }
    
    region.bedgraph <- data.frame(region.bed)
    region.bedgraph <- region.bedgraph[,c(-4,-5)]
    if(datadl) return(region.bedgraph)
    if(log2) region.bedgraph$score <- log2(region.bedgraph$score)
    
    chrom <- as.character(seqnames(region))
    chromchr <- paste(c("chr", as.character(chrom)), collapse = "")
    start <- as.integer(start(ranges(range(region))))
    end <- as.integer(end(ranges(range(region))))

    trackplot <- recordPlot()
    if(dim(region.bedgraph)[1] == 0) { # dummyplot
        plotBedpe(data.frame(), chrom, start, end, color = c("blue"), lwd = 0, plottype = "loops", heights = 0, lwdrange = c(0, 0), main = sample, adj=0)
    } else { #real plot
        plotBedgraph(region.bedgraph, chromchr, start, end, main = sample, adj=0, flip = flip)
    }
    axis(side=2,las=2,tcl=.2)
    if(showGA) labelgenome(chromchr, start, end, side = 1, scipen = 20, n = 3, scale = "Mb", line = 0.18, chromline = 0.5, scaleline = 0.5)
    return(trackplot)
}

# Primarily used for the 450k
bedgraph.trackplot <- function(file, region, ylab, sample, flip, showGA, datadl){
    region.bed <- read_delim(file, delim = " ")
    rb <- data.frame(region.bed)

    chrom <- as.character(seqnames(region))
    chromchr <- paste(c("chr", as.character(chrom)), collapse = "")
    start <- as.integer(start(ranges(range(region))))
    end <- as.integer(end(ranges(range(region))))
    
    if(datadl) return( rb[rb$seqnames==chromchr & rb$start >= start & rb$end >= end, ])

    trackplot <- recordPlot()
    plotBedgraph(rb, chromchr, start, end, main = sample, adj=0, flip = flip)
    axis(side=2,las=2,tcl=.2)
    if(showGA) labelgenome(chromchr, start, end, side = 1, scipen = 20,  n = 3, scale = "Mb", line = 0.18, chromline = 0.5, scaleline = 0.5)
    return(trackplot)
}

hicColors <- function(p) {
    if(p == 1) return(colorRampPalette(colorRamps::matlab.like2(100)))
    if(p == 2) return(colorRampPalette(c("#ffffff", "#40826D")))
    if(p == 3) return(colorRampPalette(c("#ffffff", "#9D9A96")))
    if(p == 4) return(colorRampPalette(c("#ffffff", "#2956B2")))
    if(p == 5) return(colorRampPalette(c("#ffffff", "#E34234")))
    if(p == 6) return(colorRampPalette(c("#ffffff", "#E6E6FA")))
    if(p == 7) return(colorRampPalette(c("#ffffff", "#ACE1AF")))
    if(p == 8) return(colorRampPalette(c("#ffffff", "#FF0080")))
    if(p == 9) return(colorRampPalette(c("#ffffff", "#FF9933")))
    if(p == 10) return(colorRampPalette(c("#ffffff", "#E34234")))
    if(p == 11) return(colorRampPalette(c("#ffffff", "#4B0082")))
    if(p == 12) return(colorRampPalette(c("black","blue","#1E90FF","orange","#FF8C00")))
    if(p == 13) return(colorRampPalette(c("black","blue","#1E90FF","#00BFFF","#B0E2FF")))
    if(p == 14) return(colorRampPalette(grDevices::heat.colors(100)))
    if(p == 15) return(colorRampPalette(grDevices::topo.colors(100)))
    if(p == 16) return(colorRampPalette(colorRamps::blue2red(100)))
}

hic.plot <- function(hicdata, region, sample, color, log2trans, flip, missingco, showlegend, showGA, datadl,
                     HiCmin = 0, HiCmax= 0, custMaxMin = 3, Qmin = 0, Qmax= 0){
   
    # Set up region
    chrom <- as.character(seqnames(region))
    chromchr <- paste(c("chr", as.character(chrom)), collapse = "")
    start <- as.integer(start(ranges(range(region))))
    end <- as.integer(end(ranges(range(region))))
    palette <- hicColors(color) # Hacked Sushi HiC Plot Function
    
    rows <- as.numeric(rownames(hicdata))
    cols <- as.numeric(colnames(hicdata))
    
    hicregion <- as.matrix(hicdata[which(rows >= start & rows <= end), which(cols > start & cols < end), drop=FALSE])
    if(datadl) return(data.frame(hicregion))
    
    if(log2trans) {hicregion <- log2(hicregion); hicregion[is.infinite(hicregion)] <- 0}
    
    if(dim(hicregion)[1]==0 | dim(hicregion)[2]==0){ #Nothing comes up from subsetting
        hicregion <- matrix(0)  
        colnames(hicregion) <- as.character(as.integer(start))
        rownames(hicregion) <- as.character(as.integer(end))
    }

    # determine number of bins
    rvs <- as.numeric(rownames(hicregion))
    cvs <- as.numeric(colnames(hicregion))
    min_bp <-  min(c(rvs, cvs))
    max_bp <-  max(c(rvs, cvs))
    if(length(c(diff(rvs), diff(cvs))) == 0){
        resolution <- 0
    } else {
        resolution <-  min(c(diff(rvs), diff(cvs)))
    }
    if(is.infinite(resolution)){ resolution <- max(rvs,cvs) - min(rvs,cvs)} #1x1 matrix 
    
    if(resolution != 0) {  nbins <- (max_bp-min_bp)/resolution } else { nbins <- 1 }
    
    stepsize <- abs(start - end)/(2 * nbins)
    max_z <- max(hicregion, na.rm = TRUE)
    min_z <- min(hicregion[hicregion != 0], na.rm = TRUE)    
    if(is.infinite(max_z) | is.na(max_z) | is.nan(max_z) | max_z == 0) max_z <- 10000
    if(is.infinite(min_z) | is.na(min_z) | is.nan(min_z)) min_z <- 0.0000001
    
    if(custMaxMin == 1 & (HiCmax > HiCmin)){
        max_z <- as.numeric(HiCmax)
        min_z <- as.numeric(HiCmin)
    }
    
    if(custMaxMin == 2 & (as.numeric(Qmax) > as.numeric(Qmin))){
        mreg <- melt(hicregion)
        mreg.subset <- mreg[mreg[,3] > 0  & (mreg[,1] != mreg[,2]), ]
        if(dim(mreg.subset)[1] == 0){
            max_z <- as.numeric(hicregion)
            min_z <- as.numeric(hicregion)
        } else {
            max_z <- quantile(mreg.subset[,3], as.numeric(Qmax)*0.01)
            min_z <- quantile(mreg.subset[,3], as.numeric(Qmin)*0.01)
        }
    }    
    
    # map to colors
    breaks <- seq(min_z, max_z, length.out = 100) - 0.001
    cols <- palette(length(unique(breaks)))
    
    if(missingco == "min") { cols <- c(cols[1], cols) } else { cols <- c(missingco, cols) }
    if(length(cols) == 2 ){ cols[2] <- palette(100)[100]}
    
    hicmcol <- matrix(as.character(cut(hicregion, c(-Inf, unique(breaks), Inf), labels = cols)), nrow = nrow(hicregion))
    
    # Handle flipping
    f <- 1; ylim <- c(0, 20); side <- 1
    if(flip){ f <- -1; ylim <- c(-20, 0); side <- 3}
    
    # initialize plot
    plot(1, 1, xlim = c(start, end), ylim = ylim, type = "n", xaxs = "i", yaxs = "i",
         bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = sample, adj = 0)
    if(dim(hicmcol)[1] != 1) {
    # fill plot
    h <- 20/min(40, dim(hicregion)[2]) * f
    for (rownum in (1:nrow(hicregion))) {
        y = -1*h
        x = start + (rownum * 2 * stepsize) - (stepsize * 3)
        for (colnum in (rownum:ncol(hicregion))) {
            x = x + stepsize
            y = y + h
            if((y <= 20 & f == 1) | (y >= -20 & f == -1)){
                if(colnum != rownum & y!=20*f){ # Square
                    xs = c(x - stepsize, x, x + stepsize, x, x - stepsize)
                    ys = c(y, y + h, y, y - h, y)
                } else if(y == 20 | y == -20){ #upside down triangle; at top
                    xs = c(x - stepsize, x, x + stepsize)
                    ys = c(y, y - h, y)
                } else { #basic triangle
                    xs = c(x - stepsize, x, x + stepsize)
                    ys = c(y, y + h, y)
                }
                if(rownum <= dim(hicmcol)[2] & colnum <= dim(hicmcol)[1]){
                    col <- hicmcol[colnum, rownum]
                } else {col <- cols[1]}
                polygon(xs, ys, border = NA, col = col)
            }
        }
    }
    } else {
        xs = c(start, start+stepsize, end)
        ys = c(0, f*20, 0)
        polygon(xs, ys, border = NA, col = hicmcol[1, 1])
    }
    
    if(showGA) labelgenome(chromchr, start, end, side = 1, scipen = 20, n = 3, scale = "Mb", line = 0.18, chromline = 0.5, scaleline = 0.5)
    if(min_z == max_z) min_z <- 0
    if(showlegend & !flip){
        addlegend(c(min_z, max_z), palette = palette, title="", side="right",
            bottominset=0.4, topinset=0, xoffset=-.035, labelside="left",
            width=0.025, title.offset=0.035, labels.digits=1)
    } else if(showlegend & flip) {
        addlegend(c(min_z, max_z), palette = palette, title="", side="right",
            topinset=0.4, bottominset=0.1, xoffset=-.035, labelside="left",
            width=0.025, title.offset=0.035, labels.digits=1)      
    }
}



# geneAnnotation plots the hg19/mm9 gene tracks from the cached genome loci. 
geneAnnotation <- function(y, organism, exons = FALSE, datadl) {
    chrom <- as.character(seqnames(y))
    chromchr <- paste(c("chr", as.character(chrom)), collapse = "")
    start <- as.integer(start(ranges(range(y))))
    end <- as.integer(end(ranges(range(y))))
    
    geneinfo <- data.frame()
    file <- NULL
    
    # Use cache annotation
    if(organism == "human" & exons) file <- "data/GenomeAnnotation/hg19/geneinfo-exon.rda"
    if(organism == "mouse" & exons) file <- "data/GenomeAnnotation/mm9/geneinfo-exon.rda"
    if(organism == "human" & !exons) file <- "data/GenomeAnnotation/hg19/geneinfo.rda"
    if(organism == "mouse" & !exons) file <- "data/GenomeAnnotation/mm9/geneinfo.rda"
    load(file)

    geneinfo <- geneinfo[geneinfo$chrom == chrom & geneinfo$start > start & geneinfo$stop < end,]
    if(datadl) return(geneinfo)
    
    loplot <- recordPlot()
    if(dim(geneinfo)[1] == 0){ #Dummy plot
        plotBedpe(data.frame(), chrom, start, end, color = c("blue"), lwd = 0, 
                  plottype = "loops", heights = 0, lwdrange = c(0, 0), 
                  main = "", adj=0)
    } else {
        pg <- plotGenes(geneinfo = geneinfo, chrom = chromchr, chromstart = start, 
            chromend = end, bheight = 0.1, plotgenetype = "box", 
            bentline = FALSE, labeloffset = 0.4, fontsize = 1, arrowlength = 0.025, 
            labeltext = TRUE)
    }
    #mtext(paste0("Region: ", chrom, ":", start, "-", end), outer = TRUE, line = 1)
    return(loplot)
}

