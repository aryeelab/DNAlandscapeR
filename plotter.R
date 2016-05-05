library(readr)

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

# one.loopPlot has some specialized features for plotting only 
# one sample's loops in these plots. 
one.loopPlot <- function(file, y, colorLoops = TRUE) {
    # Load the RDS
    x <- readRDS(file)
    sample <- basename(file_path_sans_ext(file))

    # Grab Regional Coordinates
    chrom <- as.character(seqnames(y))
    chromchr <- paste(c("chr", as.character(chrom)), collapse = "")
    start <- as.integer(start(ranges(range(y))))
    end <- as.integer(end(ranges(range(y))))
    
    #  Restrict the loops object to the region
    objReg <- removeSelfLoops(subsetRegion(x, y))
    
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
        lwd <- 5 * (bedPE$score/max(bedPE$score))
        
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
