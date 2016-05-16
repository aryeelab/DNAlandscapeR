masterPlotter <- function(input, dynamic.val){
    
    chia_pet_samples <- list() #Tracks samples linking with i
    chia_pet_objects <- c() #Tracks subsetted objects
    j <- 1 #Index of subsetted objects vector
    map_chia_pet.indices <- c() #maps i to order in vector of subsetted objects
    mc <- 1 #max counts
    
    #First loop initalizes the ChIA-PET max values
    for(i in input$tracks){
        i <- as.integer(i)
        if (i < 1000){ # ChIA-PET from RDS
            
            #Import object and subset
            x <- readRDS(dynamic.val$c.full[[i]])
            sample <- basename(file_path_sans_ext(dynamic.val$c.full[[i]]))
            objReg <- removeSelfLoops(subsetRegion(x, dynamic.val$region))
            
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
        if (i < 1000) {
            one.loopPlot(objReg = chia_pet_objects[[which(map_chia_pet.indices == i)]], y = dynamic.val$region,
                         sample = chia_pet_samples[[as.character(i)]], max_counts = mc)
        } else if (i < 2000) { # Track; BigWig
            bigwig.trackplot(dynamic.val$t.bw.full[[i-1000]], dynamic.val$region, "Read Depth")
        } else if (i < 3000){ # Track; Bedgraph
            bedgraph.trackplot(dynamic.val$t.bg.full[[i-2000]], dynamic.val$region, "Read Depth")
        } else if (i < 4000) { # Methyl; BigWig
            #bigwig.trackplot(dynamic.val$m.bw.full[[i-3000]], dynamic.val$region, "Methylation")
            bigwig.bumpPlot(dynamic.val$m.bw.full[[i-3000]], dynamic.val$region)
        } else if (i < 5000){ # Methyl; Bedgraph
            bedgraph.trackplot(dynamic.val$m.bg.full[[i-4000]], dynamic.val$region, "Methylation")
            #bedgraph.bumpPlot(dynamic.val$m.bw.full[[i-4000]], dynamic.val$region)
        } else {return()}
    }
    if(input$showgenes) humanAnnotation(dynamic.val$region)
    if(input$showctcf) plotCTCFregions(dynamic.val$region)
}

# one.loopPlot has some specialized features for plotting only 
# one sample's loops in these plots. The object is a loops object
# without self loops generated from the master function to determine
# the max_counts
one.loopPlot <- function(objReg, y, sample, max_counts, colorLoops = TRUE) {

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



bigwig.bumpPlot <- function(file, region, shade = TRUE){
    region.bed <- import.bw(file, which = addchr(region))
    region.bedgraph <- data.frame(region.bed)
    region.bedgraph <- region.bedgraph[,c(-4,-5)]

    chrom <- as.character(seqnames(region))
    chromchr <- paste(c("chr", as.character(chrom)), collapse = "")
    start <- as.integer(start(ranges(range(region))))
    end <- as.integer(end(ranges(range(region))))
    sample <- basename(file_path_sans_ext(file))
    
    bumpplot <- recordPlot()
    pos <- region.bedgraph$start
    y <- region.bedgraph[,4]
    cluster_id <- clusterMaker(chr=chrom, pos=pos, maxGap = 100)
    smooth <- locfitByCluster(x=pos, y=y, cluster=cluster_id, bpSpan=50)
    plot(pos, smooth$fitted, type="l", xaxt='n', ann=FALSE, main = )
    labelgenome(chromchr, start, end, side = 1, scipen = 20, 
        n = 3, scale = "Mb", line = 0.18, chromline = 0.5, scaleline = 0.5)
    mtext(sample,side=2,line=2.5,cex=1,font=2)
    if(shade) polygon(cbind(c(min(pos), pos, max(pos)), c(min(y), y, min(y))), border=NA, col="black")
    return(bumpplot)
}

bedgraph.bumpPlot <- function(file, region){
    region.bed <- read_delim(file, delim = " ")
    bedg <- data.frame(region.bed)
    region.bedgraph <- bedg[findOverlaps(GRanges(bedg), addchr(region))@from,]

    chrom <- as.character(seqnames(region))
    chromchr <- paste(c("chr", as.character(chrom)), collapse = "")
    start <- as.integer(start(ranges(range(region))))
    end <- as.integer(end(ranges(range(region))))
    sample <- basename(file_path_sans_ext(file))
    
    bumpplot <- recordPlot()
    pos <- region.bedgraph$start
    cluster_id <- clusterMaker(chr=chrom, pos=pos, maxGap = 100)
    smooth <- locfitByCluster(x=pos, y=region.bedgraph[,4], cluster=cluster_id, bpSpan=50)
    plot(pos, smooth$fitted, type="l", xaxt='n', ann=FALSE)
    labelgenome(chromchr, start, end, side = 1, scipen = 20, 
        n = 3, scale = "Mb", line = 0.18, chromline = 0.5, scaleline = 0.5)
    mtext("Methylation",side=2,line=2.5,cex=1,font=2)
    points(pos, region.bedgraph[,4])
    return(bumpplot)
}

# plots bigwig data for specified file/region; annotates 'ylab' 
bigwig.trackplot <- function(file, region, ylab){
    region.bed <- import.bw(file, which = addchr(region))
    region.bedgraph <- data.frame(region.bed)
    region.bedgraph <- region.bedgraph[,c(-4,-5)]
    
    chrom <- as.character(seqnames(region))
    chromchr <- paste(c("chr", as.character(chrom)), collapse = "")
    start <- as.integer(start(ranges(range(region))))
    end <- as.integer(end(ranges(range(region))))
    sample <- basename(file_path_sans_ext(file))
    
    trackplot <- recordPlot()
    plotBedgraph(region.bedgraph, chromchr, start, end, 
                 main = sample, adj=0)
    mtext(ylab,side=2,line=2.5,cex=1,font=2)
    axis(side=2,las=2,tcl=.2)
    labelgenome(chromchr, start, end, side = 1, scipen = 20, 
                n = 3, scale = "Mb", line = 0.18, chromline = 0.5, scaleline = 0.5)
    return(trackplot)
}

# plots bedgraph data for specified file/region; annotates 'ylab'
bedgraph.trackplot <- function(file, region, ylab){
    region.bed <- read_delim(file, delim = " ")
    region.bedgraph <- data.frame(region.bed)

    chrom <- as.character(seqnames(region))
    chromchr <- paste(c("chr", as.character(chrom)), collapse = "")
    start <- as.integer(start(ranges(range(region))))
    end <- as.integer(end(ranges(range(region))))
    sample <- basename(file_path_sans_ext(file))
    
    trackplot <- recordPlot()
    plotBedgraph(region.bedgraph, chromchr, start, end, 
                 main = sample, adj=0)
    mtext(ylab,side=2,line=2.5,cex=1,font=2)
    axis(side=2,las=2,tcl=.2)
    labelgenome(chromchr, start, end, side = 1, scipen = 20, 
                n = 3, scale = "Mb", line = 0.18, chromline = 0.5, scaleline = 0.5)
    return(trackplot)
}

# humanAnnotation plots the human gene tracks from the cached genome loci. 
humanAnnotation <- function(y) {
    chrom <- as.character(seqnames(y))
    chromchr <- paste(c("chr", as.character(chrom)), collapse = "")
    start <- as.integer(start(ranges(range(y))))
    end <- as.integer(end(ranges(range(y))))
    
    # Use cache annotation
    rda <- paste(system.file("rda", package = "diffloop"), "geneinfo.h.rda", sep = "/")
    load(rda)
    geneinfo <- geneinfo[geneinfo$chrom == chrom & geneinfo$start > start - 10000 & geneinfo$stop < end + 10000,]

    loplot <- recordPlot()
    pg = plotGenes(geneinfo = geneinfo, chrom = chromchr, chromstart = start, 
        chromend = end, bheight = 0.1, plotgenetype = "box", 
        bentline = FALSE, labeloffset = 0.4, fontsize = 1, arrowlength = 0.025, 
        labeltext = TRUE)
    mtext(paste0("Region: ", chrom, ":", start, "-", end), outer = TRUE, 
        line = 1)
    return(loplot)
}


plotCTCFregions <- function(y) {
    chrom <- as.character(seqnames(y))
    chromchr <- paste(c("chr", as.character(chrom)), collapse = "")
    start <- as.integer(start(ranges(range(y))))
    end <- as.integer(end(ranges(range(y))))
    
    ctcf <- data.frame(read_delim("data/CTCF-regions.bed", delim = "\t"))
    ctcf.small <- ctcf[findOverlaps(GRanges(ctcf), addchr(y))@from,]
    
    loplot <- recordPlot()
    pg = plotBed(beddata = ctcf.small, chrom = chromchr, chromstart = start, 
        chromend = end, labeltext = TRUE)
    labelgenome(chromchr, start, end, side = 1, scipen = 20, 
                    n = 3, scale = "Mb", line = 0.18, chromline = 0.5, scaleline = 0.5)
    mtext(paste0("Region: ", chrom, ":", start, "-", end), outer = TRUE, 
        line = 1)
    mtext("CTCF Sites",side=2,line=2.5,cex=1,font=2)

    return(loplot)
}
