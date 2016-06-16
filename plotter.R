# Main functions for plotting the various tracks. A brief overview follows.

# masterPlotter is called from the server.R script and performs these operations:
# 1) subsets all ChIA-PET objects to determine the max counts for normalization
# 2) plots tracks calling later functions based on the indices.

masterPlotter <- function(input, dynamic.val){
    
    chia_pet_samples <- list() #Tracks samples linking with i
    chia_pet_objects <- c() #Tracks subsetted objects
    j <- 1 #Index of subsetted objects vector
    map_chia_pet.indices <- c() #maps i to order in vector of subsetted objects
    mc <- 1 #max counts
    one_anchor_samples <- list() #Tracks samples linking with i
    
    #First loop initalizes the ChIA-PET max values
    for(i in input$tracks){
        i <- as.integer(i)
        if (i < 1000000){ # ChIA-PET from RDS
            
            #Import object and subset
            file.conn <- dynamic.val$c.full[[i]]
            if(grepl("amazonaws", file.conn)){ x <- readRDS(gzcon(url(file.conn)))
            } else { x <- readRDS(file.conn) }
            sample <- names(dynamic.val$c.list)[i]
            objReg <- removeSelfLoops(.subsetRegion.quick(x, dynamic.val$region, nanchors = 2))
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
    
    #Second loop does all the plotting
    for(i in input$tracks){ 
        i <- as.integer(i)
        if (i < 1000000) {
            oa <- try(one_anchor_samples[[which(map_chia_pet.indices == i)]], silent = TRUE)
            if("try-error" %in% class(oa)) oa <- NULL
            one.loopPlot(objReg = chia_pet_objects[[which(map_chia_pet.indices == i)]], y = dynamic.val$region,
                         sample = chia_pet_samples[[as.character(i)]], max_counts = mc, oneAnchor = oa)
        } else if (i < 2000000) { # Track; BigWig
            t <- i - 1000000
            sample <- names(dynamic.val$t.bw.list)[t]
            bigwig.trackplot(dynamic.val$t.bw.full[[t]], dynamic.val$region, input$smoother, FUN = input$FUN, "Depth", sample = sample, log2 = input$log2BW)
        } else if (i < 3000000){ # Track; Bedgraph
            t <- i - 2000000
            sample <- names(dynamic.val$t.bg.list)[t]
            bedgraph.trackplot(dynamic.val$t.bg.full[[t]], dynamic.val$region, "Depth", sample = sample)
        } else if (i < 4000000) { # Methyl; BigWig
            t <- i - 3000000
            sample <- names(dynamic.val$m.bw.list)[t]
            bigwig.bumpPlot(dynamic.val$m.bw.full[[t]], dynamic.val$region, sample = sample)
        } else if (i < 5000000){ # Methyl; Bedgraph
            t <- i - 4000000
            sample <- names(dynamic.val$m.bg.list)[t]
            bedgraph.trackplot(dynamic.val$m.bg.full[[t]], dynamic.val$region, "Methylation", sample = sample)
        } else if (i < 6000000){ # Hi-C Plot for Stuff that's preloaded on the server
            t <- i - 5000000
            sample.hic <- names(dynamic.val$i.list)[t]
            sample <- gsub("-HiC", "", sample.hic)
            fs <- g_h.i.full
            res <- as.character(input[[paste0(sample, "HiCRes")]])
            chrom <- paste0("chr", as.character(seqnames(dynamic.val$region)))
            file <- fs[grepl(paste0(chrom, ".rds"), fs) & grepl(res, fs) & grepl(sample, fs)]
            hicdata <- readRDS(gzcon(url(file)))
            hic.plot(hicdata, dynamic.val$region, sample = sample.hic, color = input$HiCcolor, log2trans = input$log2hic)
        } else if (i < 7000000) { # Local Hi-C plot    
            t <- i - 6000000
            chrom <- paste0("chr", as.character(seqnames(dynamic.val$region)))
            list.dat <- readRDS(dynamic.val$i.l.full[t])
            hicdata <- list.dat[[chrom]]
            sample <- names(dynamic.val$i.l.list)[t]
            hic.plot(hicdata, dynamic.val$region, sample = sample, color = input$HiCcolor, log2trans = input$log2hic)
        } else {return()}
    }
    e <- ifelse(input$showgenes == 2, TRUE, FALSE)
    if(input$showgenes > 0 & input$organism == 1) geneAnnotation(dynamic.val$region, "human", exons = e)
    if(input$showgenes > 0 & input$organism == 2) geneAnnotation(dynamic.val$region, "mouse", exons = e)
}

# one.loopPlot has some specialized features for plotting only 
# one sample's loops in these plots. The object is a loops object
# without self loops generated from the master function to determine
# the max_counts
one.loopPlot <- function(objReg, y, sample, max_counts, colorLoops = TRUE, oneAnchor = NULL) {

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
        bedPE <- data.frame(LA, RA, name, score, strand_1, strand_2, sample)

        w <- loopWidth(objReg)
        h <- sqrt(w/max(w))
        lwd <- 5 * (bedPE$score/max_counts)
        
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
        plotBedpe(bedPE, chrom, start, end, color = cs, lwd = lwd, 
                  plottype = "loops", heights = h, lwdrange = c(0, 5), 
                  main = sample, adj=0)
        labelgenome(chromchr, start, end, side = 1, scipen = 20, 
                    n = 3, scale = "Mb", line = 0.18, chromline = 0.5, scaleline = 0.5)
        return(loplot)
    } else {
        # Return dummy plot
        loplot <- recordPlot()
        plotBedpe(data.frame(), chrom, start, end, color = c("blue"), lwd = 0, 
                  plottype = "loops", heights = 0, lwdrange = c(0, 0), 
                  main = sample, adj=0)
        labelgenome(chromchr, start, end, side = 1, scipen = 20, 
                    n = 3, scale = "Mb", line = 0.18, chromline = 0.5, scaleline = 0.5)
        return(loplot)   
    }
}

# bigwig.bumpPlot is used for methylation
bigwig.bumpPlot <- function(file, region, smoother = 0, shade = TRUE, sample){
    region.bed <- import.bw(file, which = addchr(region))
    
    region.bedgraph <- data.frame(region.bed)
    region.bedgraph <- region.bedgraph[,c(-4,-5)]

    chrom <- as.character(seqnames(region))
    chromchr <- paste(c("chr", as.character(chrom)), collapse = "")
    start <- as.integer(start(ranges(range(region))))
    end <- as.integer(end(ranges(range(region))))

    bumpplot <- recordPlot()
    pos <- region.bedgraph$start
    y <- region.bedgraph[,4]
    cluster_id <- clusterMaker(chr=chrom, pos=pos, maxGap = 100)
    smooth <- locfitByCluster(x=pos, y=y, cluster=cluster_id, bpSpan = smoother + 50)
    plot(pos, smooth$fitted, type="l", xaxt='n',bty = "n",xaxs="i",yaxs="i",main=sample,adj=0,ylab="")
    labelgenome(chromchr, start, end, side = 1, scipen = 20, 
        n = 3, scale = "Mb", line = 0.18, chromline = 0.5, scaleline = 0.5)
    # mtext(sample,side=2,line=2.5,cex=1,font=2)
    if(shade) polygon(cbind(c(min(pos), pos, max(pos)), c(min(y), y, min(y))), border=NA, col="black")
    return(bumpplot)
}

# bigwig.trackplot is used for most epigenetic peaks
bigwig.trackplot <- function(file, region, smoother, FUN, ylab, sample, log2){
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
    if(log2) region.bedgraph$score <- log2(region.bedgraph$score)
    
    chrom <- as.character(seqnames(region))
    chromchr <- paste(c("chr", as.character(chrom)), collapse = "")
    start <- as.integer(start(ranges(range(region))))
    end <- as.integer(end(ranges(range(region))))

    trackplot <- recordPlot()
    plotBedgraph(region.bedgraph, chromchr, start, end, 
                 main = sample, adj=0)
    #mtext(ylab,side=2,line=2.5,cex=1,font=2)
    axis(side=2,las=2,tcl=.2)
    labelgenome(chromchr, start, end, side = 1, scipen = 20, 
                n = 3, scale = "Mb", line = 0.18, chromline = 0.5, scaleline = 0.5)
    return(trackplot)
}

# Primarily used for the 450k
bedgraph.trackplot <- function(file, region, ylab, sample){
    region.bed <- read_delim(file, delim = " ")
    region.bedgraph <- data.frame(region.bed)

    chrom <- as.character(seqnames(region))
    chromchr <- paste(c("chr", as.character(chrom)), collapse = "")
    start <- as.integer(start(ranges(range(region))))
    end <- as.integer(end(ranges(range(region))))

    trackplot <- recordPlot()
    plotBedgraph(region.bedgraph, chromchr, start, end, main = sample, adj=0)
    axis(side=2,las=2,tcl=.2)
    labelgenome(chromchr, start, end, side = 1, scipen = 20, 
                n = 3, scale = "Mb", line = 0.18, chromline = 0.5, scaleline = 0.5)
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

hic.plot <- function(hicdata, region, sample, color, log2trans){
    chrom <- as.character(seqnames(region))
    chromchr <- paste(c("chr", as.character(chrom)), collapse = "")
    start <- as.integer(start(ranges(range(region))))
    end <- as.integer(end(ranges(range(region))))
    
    # Hacked Sushi HiC Plot Function
    palette <- hicColors(color)
    rows <- as.numeric(rownames(hicdata))
    cols <- as.numeric(colnames(hicdata))
    
    hicregion <- as.matrix(hicdata[which(rows >= start & rows <= end), which(cols >= start & cols <= end)])
    if(log2trans) {
        hicregion <- log2(hicregion)
        hicregion[hicregion < 0] <- 0
    }
    
    # determine number of bins
    nbins <- nrow(hicregion)
    stepsize <- abs(start - end)/(2 * nbins)
    max_z <- max(hicregion, na.rm = TRUE)
    min_z <- min(hicregion[hicregion > 0], na.rm = TRUE)    
        
    # map to colors
    breaks <- seq(min_z, max_z, length.out = 100)
    cols <- palette(length(breaks) + 1)
    if(length(unique(breaks)) == 1){
        hicmcol <- matrix(max(cols), nrow = 1, ncol = 1)
    } else {
        hicmcol <- matrix(as.character(cut(hicregion, c(-Inf, breaks, Inf), labels = cols)), nrow = nrow(hicregion))
    }
    # initialize plot
    plot(1, 1, xlim = c(start, end), ylim = c(0, 20), type = "n", xaxs = "i", yaxs = "i",
         bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = sample, adj = 0)

    # fill plot
    h <- 20/min(40, dim(hicregion)[2])
    for (rownum in (1:nrow(hicregion))) {
        y = -1*h
        x = start + (rownum * 2 * stepsize) - (stepsize * 2)
        for (colnum in (rownum:ncol(hicregion))) {
            x = x + stepsize
            y = y + h
            if(y <= 20){
                if(colnum != rownum & y!=20){ # Square
                    xs = c(x - stepsize, x, x + stepsize, x, x - stepsize)
                    ys = c(y, y + h, y, y - h, y)
                } else if(y == 20){ #upside down triangle
                    xs = c(x - stepsize, x, x + stepsize)
                    ys = c(y, y - h, y)
                } else {
                    xs = c(x - stepsize, x, x + stepsize)
                    ys = c(y, y + h, y)
                }
                polygon(xs, ys, border = NA, col = hicmcol[colnum, rownum])
            }
        }
    }
    labelgenome(chromchr, start, end, n=4, scale="Mb",edgeblankfraction=0.20)
    if(min_z == max_z) min_z <- 0
    addlegend(c(min_z, max_z), palette = palette, title="", side="right",
        bottominset=0.4, topinset=0, xoffset=-.035, labelside="left",
        width=0.025, title.offset=0.035)
}

# geneAnnotation plots the hg19/mm9 gene tracks from the cached genome loci. 
geneAnnotation <- function(y, organism, exons = FALSE) {
    chrom <- as.character(seqnames(y))
    chromchr <- paste(c("chr", as.character(chrom)), collapse = "")
    start <- as.integer(start(ranges(range(y))))
    end <- as.integer(end(ranges(range(y))))
    
    geneinfo <- data.frame()
    
    # Use cache annotation
    if(organism == "human" & exons) load("data/GenomeAnnotation/hg19/geneinfo-exon.rda")
    if(organism == "mouse" & exons) load("data/GenomeAnnotation/mm9/geneinfo-exon.rda")
    if(organism == "human" & !exons) load("data/GenomeAnnotation/hg19/geneinfo.rda")
    if(organism == "mouse" & !exons) load("data/GenomeAnnotation/mm9/geneinfo.rda")
    
    geneinfo <- geneinfo[geneinfo$chrom == chrom & geneinfo$start > start & geneinfo$stop < end,]

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

