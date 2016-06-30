# Function to import another bucket #

importAmazonAWSBucket <- function(newBucket, dynamic.val){
    
    # Initialize data from Amazon
    amazon <- paste0("http://s3.amazonaws.com/", newBucket)
    t <- unlist(get_bucket(bucket = newBucket))
    amazon.filenames <- paste(amazon, t[grep("data", t)], sep = "/")

    # Append Amazon data
    g_h.c.full <-  c(dynamic.val$h.c.full, amazon.filenames[grepl("data/human/loops/.{1,}", amazon.filenames)])
    g_h.t.files <- c(dynamic.val$h.t.files, amazon.filenames[grepl("data/human/tracks/.{1,}", amazon.filenames)])
    g_h.m.files <- c(dynamic.val$h.m.files, amazon.filenames[grepl("data/human/methylation/.{1,}", amazon.filenames)])
    
    # Append Amazon HiC data
    i.temp <- amazon.filenames[grepl("data/human/hic/.{1,}", amazon.filenames)]
    i.base <- basename(i.temp)
    amazon.hic.samples <- basename(i.base[!grepl(".rds", i.base) & !grepl("000", i.base)])
    g_h.i.samples <- c(gsub("-HiC", "", names(dynamic.val$h.i.list)), amazon.hic.samples)
    res.temp2 <-  unlist(strsplit(basename(i.temp[grepl(".rds", i.temp) & grepl("_", i.temp)]), "-chr"))
    res.temp <- unique(res.temp2[!grepl(".rds", res.temp2) ])
    g_h.i.res <- c(dynamic.val$h.i.res, lapply(amazon.hic.samples, function(t){
        opts <- unlist(strsplit(res.temp[grepl(t, res.temp)],split="_"))
        unique(opts[grep("000", opts)])
    }))
    names(g_h.i.res) <- g_h.i.samples
    g_h.i.full <- c(dynamic.val$h.i.full, i.temp[grepl(".rds", i.temp)])

    # From 1-1,000,000-- ChIA-PET loops objects
    if(length(g_h.c.full) != 0){
        g_h.c.names <- basename(file_path_sans_ext(g_h.c.full))
        g_h.c.list <- as.list(seq(1, length(g_h.c.names), by = 1) + 0)
        names(g_h.c.list) <- g_h.c.names
    } else { g_h.c.list <- list(); g_h.c.full <- list()}
    
    bigwig <- c(".bw", ".bigwig", ".bigWig")
    # From 1,000,001-2,000,000-- ReadDepth Tracks-- bigwig
    if(length(g_h.t.files) != 0){
        g_h.t.bw.full <- g_h.t.files[as.logical(rowSums(sapply(bigwig, grepl, g_h.t.files)))]
        g_h.t.bw.names <- basename(file_path_sans_ext(g_h.t.bw.full))
        g_h.t.bw.list <- as.list(seq(1, length(g_h.t.bw.names), by = 1) + 1000000)
        names(g_h.t.bw.list) <- g_h.t.bw.names
    } else { g_h.t.bw.list <- list(); g_h.t.bw.full <- list()}
    
    # From 2,000,001-3,000,000-- ReadDepth Tracks-- Bedgraph
    g_h.t.bg.full <- g_h.t.files[grep(".bedgraph", g_h.t.files, fixed=T)]
    if(length(g_h.t.bg.full) != 0){
        g_h.t.bg.names <- basename(file_path_sans_ext(g_h.t.bg.full))
        g_h.t.bg.list <- as.list(seq(1, length(g_h.t.bg.names), by = 1) + 2000000)
        names(g_h.t.bg.list) <- g_h.t.bg.names
    } else { g_h.t.bg.list <- list(); g_h.t.bg.full <- list() }
    
    # From 3,000,001-4,000,000-- Methylation Tracks-- bigwig
    if(length(g_h.m.files) != 0){
        g_h.m.bw.full <- g_h.m.files[as.logical(rowSums(sapply(bigwig, grepl, g_h.m.files)))]
        g_h.m.bw.names <- basename(file_path_sans_ext(g_h.m.bw.full))
        g_h.m.bw.list <- as.list(seq(1, length(g_h.m.bw.names), by = 1) + 3000000)
        names(g_h.m.bw.list) <- g_h.m.bw.names
    } else { g_h.m.bw.full <- list(); g_h.m.bw.list <- list(); g_h.m.bw.full <- list()}
    
    # From 4,000,001-5,000,000-- Methylation Tracks-- Bedgraph
    g_h.m.bg.full <- g_h.m.files[grep(".bedgraph", g_h.m.files, fixed=T)]
    if(length(g_h.m.bg.full) != 0){
        g_h.m.bg.names <- basename(file_path_sans_ext(g_h.m.bg.full))
        g_h.m.bg.list <- as.list(seq(1, length(g_h.m.bg.names), by = 1) + 4000000)
        names(g_h.m.bg.list) <- g_h.m.bg.names
    } else { g_h.m.bg.list <- list(); g_h.m.bg.full <- list()}
    
    # From 5,000,001-6,000,000-- HiC Tracks-- .rds
    if(length(g_h.i.samples) != 0){
        g_h.i.list <- as.list(seq(1, length(g_h.i.samples), by = 1) + 5000000)
        names(g_h.i.list) <- paste(g_h.i.samples, "-HiC", sep="")
    } else {g_h.i.list <- list()}
    
    # From 6,000,001-7,000,000-- Local HiC Tracks-- .rds
    # Do Nothing; just keeping track. 
    
    g_h.f.list <- append(g_h.c.list, append(append(g_h.t.bw.list, g_h.t.bg.list), append(append(g_h.m.bw.list, g_h.m.bg.list), g_h.i.list)))
    
    # Now order it
    h.d <- data.frame(unlist(g_h.f.list))
    h.d$names <- rownames(h.d)
    h.d <- h.d[order(rownames(h.d)), ]
    g_h.f.list <- list()
    for(k in 1:dim(h.d)[1]){
        x <- list(h.d[k,1])
        names(x) <- (h.d[k,2])
        g_h.f.list <- append(g_h.f.list, x)
    }
    
    # Append amazon data
    g_m.c.full <- c(g_m.c.full, amazon.filenames[grepl("data/mouse/loops/.{1,}", amazon.filenames)])
    g_m.t.files <- c(g_m.t.files, amazon.filenames[grepl("data/mouse/tracks/.{1,}", amazon.filenames)])
    g_m.m.files <- c(g_m.m.files, amazon.filenames[grepl("data/mouse/methylation/.{1,}", amazon.filenames)])
    
    # Append Amazon HiC data
    i.temp <- amazon.filenames[grepl("data/mouse/hic/.{1,}", amazon.filenames)]
    i.base <- basename(i.temp)
    amazon.hic.samples <- basename(i.base[!grepl(".rds", i.base) & !grepl("000", i.base)])
    g_m.i.samples <- c(gsub("-HiC", "", names(dynamic.val$m.i.list)), amazon.hic.samples)
    res.temp2 <-  unlist(strsplit(basename(i.temp[grepl(".rds", i.temp) & grepl("_", i.temp)]), "-chr"))
    res.temp <- unique(res.temp2[!grepl(".rds", res.temp2) ])
    g_m.i.res <- c(dynamic.val$m.i.res, lapply(amazon.hic.samples, function(t){
        opts <- unlist(strsplit(res.temp[grepl(t, res.temp)],split="_"))
        unique(opts[grep("000", opts)])
    }))
    names(g_m.i.res) <- g_m.i.samples
    g_m.i.full <- c(dynamic.val$m.i.full, i.temp[grepl(".rds", i.temp)])
    
    # From 1-1,000,000-- ChIA-PET loops objects
    if(length(g_m.c.full) != 0){
        g_m.c.names <- basename(file_path_sans_ext(g_m.c.full))
        g_m.c.list <- as.list(seq(1, length(g_m.c.names), by = 1) + 0)
        names(g_m.c.list) <- g_m.c.names
    } else { g_m.c.list <- list(); g_m.c.full <- list()}
    
    bigwig <- c(".bw", ".bigwig", ".bigWig")
    # From 1,000,001-2,000,000-- ReadDepth Tracks-- bigwig
    if(length(g_m.t.files) != 0){
        g_m.t.bw.full <- g_m.t.files[as.logical(rowSums(sapply(bigwig, grepl, g_m.t.files)))]
        g_m.t.bw.names <- basename(file_path_sans_ext(g_m.t.bw.full))
        g_m.t.bw.list <- as.list(seq(1, length(g_m.t.bw.names), by = 1) + 1000000)
        names(g_m.t.bw.list) <- g_m.t.bw.names
    } else {g_m.t.bw.full <- c();  g_m.t.bw.list <- list(); g_m.t.bw.full <- list()}
    
    # From 2,000,001-3,000,000-- ReadDepth Tracks-- Bedgraph
    g_m.t.bg.full <- g_m.t.files[grep(".bedgraph", g_m.t.files, fixed=T)]
    if(length(g_m.t.bg.full) != 0){
        g_m.t.bg.names <- basename(file_path_sans_ext(g_m.t.bg.full))
        g_m.t.bg.list <- as.list(seq(1, length(g_m.t.bg.names), by = 1) + 2000000)
        names(g_m.t.bg.list) <- g_m.t.bg.names
    } else { g_m.t.bg.list <- list(); g_m.t.bg.full <- list() }
    
    # From 3,000,001-4,000,000-- Methylation Tracks-- bigwig
    if(length(g_m.m.files) != 0){
        temp <- matrix(sapply(bigwig, grepl, g_m.m.files), ncol = 3)
        g_m.m.bw.full <- g_m.m.files[as.logical(rowSums(temp))]
        g_m.m.bw.names <- basename(file_path_sans_ext(g_m.m.bw.full))
        g_m.m.bw.list <- as.list(seq(1, length(g_m.m.bw.names), by = 1) + 3000000)
        names(g_m.m.bw.list) <- g_m.m.bw.names
    } else { g_m.m.bw.full <- list(); g_m.m.bw.list <- list(); g_m.m.bw.full <- list()}
    
    # From 4,000,001-5,000,000-- Methylation Tracks-- Bedgraph
    g_m.m.bg.full <- g_m.m.files[grep(".bedgraph", g_m.m.files, fixed=T)]
    if(length(g_m.m.bg.full) != 0){
        g_m.m.bg.names <- basename(file_path_sans_ext(g_m.m.bg.full))
        g_m.m.bg.list <- as.list(seq(1, length(g_m.m.bg.names), by = 1) + 4000000)
        names(g_m.m.bg.list) <- g_m.m.bg.names
    } else { g_m.m.bg.list <- list(); g_m.m.bg.full <- list()}
    
    # From 5,000,001-6,000,000-- HiC Tracks-- .rds -- not currently supported
    if(length(g_m.i.samples) != 0){
        g_m.i.list <- as.list(seq(1, length(g_m.i.samples), by = 1) + 5000000)
        names(g_m.i.list) <- paste(g_m.i.samples, "-HiC", sep="")
    } else {g_m.i.list <- list()}
    
    # From 6,000,001-7,000,000-- Local HiC Tracks-- .rds
    # Do Nothing; just keeping track. 
    
    g_m.f.list <- append(g_m.c.list, append(append(g_m.t.bw.list, g_m.t.bg.list), append(append(g_m.m.bw.list, g_m.m.bg.list), g_m.i.list)))
    
    # Now order it
    m.d <- data.frame(unlist(g_m.f.list))
    m.d$names <- rownames(m.d)
    m.d <- m.d[order(rownames(m.d)), ]
    m.f.list <- list()
    for(k in 1:dim(m.d)[1]){
        x <- list(m.d[k,1])
        names(x) <- (m.d[k,2])
        g_m.f.list <- append(g_m.f.list, x)
    }
    
    # Update dynamic values
    
    # Update human variables
	dynamic.val$h.f.list    = g_h.f.list
    dynamic.val$h.c.full    = g_h.c.full
    dynamic.val$h.t.files   = g_h.t.files
    dynamic.val$h.m.files   = g_h.m.files
    dynamic.val$h.t.bw.full = g_h.t.bw.full
    dynamic.val$h.t.bg.full = g_h.t.bg.full
    dynamic.val$h.m.bw.full = g_h.m.bw.full
    dynamic.val$h.m.bg.full = g_h.m.bg.full
    dynamic.val$h.i.full    = g_h.i.full
    dynamic.val$h.i.res     = g_h.i.res
    dynamic.val$h.i.l.full  = NULL
    dynamic.val$h.c.list    = g_h.c.list
    dynamic.val$h.t.bw.list = g_h.t.bw.list
    dynamic.val$h.t.bg.list = g_h.t.bg.list
    dynamic.val$h.m.bw.list = g_h.m.bw.list
    dynamic.val$h.m.bg.list = g_h.m.bg.list
    dynamic.val$h.i.list    = g_h.i.list
    dynamic.val$h.i.l.list  = NULL
    
    # Update mouse variables
    dynamic.val$m.f.list    = g_m.f.list
    dynamic.val$m.t.files   = g_m.t.files
    dynamic.val$m.m.files   = g_m.m.files
    dynamic.val$m.c.full    = g_m.c.full
    dynamic.val$m.t.bw.full = g_m.t.bw.full
    dynamic.val$m.t.bg.full = g_m.t.bg.full
    dynamic.val$m.m.bw.full = g_m.m.bw.full
    dynamic.val$m.m.bg.full = g_m.m.bg.full
    dynamic.val$m.i.full    = g_m.i.full
    dynamic.val$m.i.res     = g_m.i.res
    dynamic.val$m.i.l.full  = NULL
    dynamic.val$m.c.list    = g_m.c.list
    dynamic.val$m.t.bw.list = g_m.t.bw.list
    dynamic.val$m.t.bg.list = g_m.t.bg.list
    dynamic.val$m.m.bw.list = g_m.m.bw.list
    dynamic.val$m.m.bg.list = g_m.m.bg.list
    dynamic.val$m.i.list    = g_m.i.list
    dynamic.val$m.i.l.list  = NULL
    return(dynamic.val)
}
