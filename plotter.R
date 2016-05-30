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
    
    #First loop initalizes the ChIA-PET max values
    for(i in input$tracks){
        i <- as.integer(i)
        if (i < 1000000){ # ChIA-PET from RDS
            
            #Import object and subset
            x <- readRDS(dynamic.val$c.full[[i]])
            sample <- names(dynamic.val$c.list)[i]
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
        if (i < 1000000) {
            one.loopPlot(objReg = chia_pet_objects[[which(map_chia_pet.indices == i)]], y = dynamic.val$region,
                         sample = chia_pet_samples[[as.character(i)]], max_counts = mc)
        } else if (i < 2000000) { # Track; BigWig
            t <- i - 1000000
            sample <- names(dynamic.val$t.bw.list)[t]
            bigwig.trackplot(dynamic.val$t.bw.full[[t]], dynamic.val$region, "Depth", sample = sample)
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
            sample <- names(dynamic.val$t.bg.list)[t]
            bedgraph.trackplot(dynamic.val$m.bg.full[[t]], dynamic.val$region, "Methylation", sample = sample)
        } else {return()}
    }
    
    if(input$showgenes & input$organism == 1) geneAnnotation(dynamic.val$region, "human", input$plotGenes)
    if(input$showgenes & input$organism == 2) geneAnnotation(dynamic.val$region, "mouse", input$plotGenes)
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

# bigwig.bumpPlot is used for methylation
bigwig.bumpPlot <- function(file, region, shade = TRUE, sample){
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
    smooth <- locfitByCluster(x=pos, y=y, cluster=cluster_id, bpSpan=50)
    plot(pos, smooth$fitted, type="l", xaxt='n', ann=FALSE, bty = "n",xaxs="i",yaxs="i")
    labelgenome(chromchr, start, end, side = 1, scipen = 20, 
        n = 3, scale = "Mb", line = 0.18, chromline = 0.5, scaleline = 0.5)
    # mtext(sample,side=2,line=2.5,cex=1,font=2)
    if(shade) polygon(cbind(c(min(pos), pos, max(pos)), c(min(y), y, min(y))), border=NA, col="black")
    return(bumpplot)
}

# bigwig.trackplot is used for most epigenetic peaks
bigwig.trackplot <- function(file, region, ylab, sample){
    region.bed <- import.bw(file, which = addchr(region))
    region.bedgraph <- data.frame(region.bed)
    region.bedgraph <- region.bedgraph[,c(-4,-5)]
    
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
    plotBedgraph(region.bedgraph, chromchr, start, end, 
                 main = sample, adj=0)
    mtext(ylab,side=2,line=2.5,cex=1,font=2)
    axis(side=2,las=2,tcl=.2)
    labelgenome(chromchr, start, end, side = 1, scipen = 20, 
                n = 3, scale = "Mb", line = 0.18, chromline = 0.5, scaleline = 0.5)
    return(trackplot)
}

# geneAnnotation plots the hg19/mm9 gene tracks from the cached genome loci. 
geneAnnotation <- function(y, organism, plotGenes) {
    chrom <- as.character(seqnames(y))
    chromchr <- paste(c("chr", as.character(chrom)), collapse = "")
    start <- as.integer(start(ranges(range(y))))
    end <- as.integer(end(ranges(range(y))))
    
    geneinfo <- data.frame()
    
    # Use cache annotation
    if(organism == "human") load("data/GenomeAnnotation/hg19/geneinfo.rda")
    if(organism == "mouse") load("data/GenomeAnnotation/mm9/geneinfo.rda")
    
    geneinfo <- geneinfo[geneinfo$chrom == chrom & geneinfo$start > start & geneinfo$stop < end,]
    geneinfo <- geneinfo[geneinfo$gene %in% plotGenes, ]

    loplot <- recordPlot()
    if(dim(geneinfo)[1] == 0){
        plotBedpe(data.frame(), chrom, start, end, color = c("blue"), lwd = 0, 
                  plottype = "loops", heights = 0, lwdrange = c(0, 0), 
                  main = "", adj=0)
    } else {
        pg <- plotGenes(geneinfo = geneinfo, chrom = chromchr, chromstart = start, 
            chromend = end, bheight = 0.1, plotgenetype = "box", 
            bentline = FALSE, labeloffset = 0.4, fontsize = 1, arrowlength = 0.025, 
            labeltext = TRUE)
    }
    mtext(paste0("Region: ", chrom, ":", start, "-", end), outer = TRUE, 
        line = 1)
    return(loplot)
}

