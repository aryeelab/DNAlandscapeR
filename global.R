# Global variables for both server.R and ui.R to reference
library(tools)
library(shiny)
library(shinyBS)
library(shinythemes)
library(ggplot2)
library(GenomicRanges)
library(diffloop)
library(Sushi)
library(foreach)
library(rtracklayer)
library(DT)
library(grid)
library(readr)
library(bumphunter)
library(shinyFiles)
library(markdown)
library(knitr)
library(gsubfn)    
library(RCurl)
library(Matrix)

default_chr <- "9"
default_start <- 21912689
default_end <- 22216233

uploadchoices <- list("Loops/.rds" = 1,
                      "ReadDepth/.bigWig" = 2,
                      "ReadDepth/.bedgraph" = 3,
                      "Methylation/.bigWig" = 4,
                      "Methylation/.bedgraph" = 5,
                      "HiC (Individual Chromosome)/.rds" = 6)

# Variable names are coded as follows
# g_ is a global variable
# h. is human; m. is mouse
# c. is ChIA-PET data
# i. is HiC data
# t. and m. (in the third position) are track and methylation files; f. is full (all)
# bw. and bw. are bigwig and bedgraph files
# full is the path to the file


# Initialize data from Amazon
amazon <- "http://s3.amazonaws.com/dnalandscaper"
xmlDat <- getURL(amazon, ftp.use.epsv = FALSE, dirlistonly = TRUE)
amazon.filenames <- gsubfn::strapplyc(xmlDat, "<Contents><Key>(.*?)</Key>", simplify = c)
amazon.filenames <- paste(amazon, amazon.filenames, sep = "/")


## HUMAN INITIALIZATION ##

# Import locally hosted data file names
g_h.c.full <- list.files("data/human/loops", full.names = TRUE)
g_h.t.files <- list.files("data/human/tracks", full.names = TRUE)
g_h.m.files <- list.files("data/human/methylation", full.names = TRUE)
g_h.i.full <- list.files("data/human/hic", full.names = TRUE)

# Append amazon data
g_h.c.full <- c(g_h.t.files, amazon.filenames[grepl("data/human/loops/.{1,}", amazon.filenames)])
g_h.t.files <- c(g_h.t.files, amazon.filenames[grepl("data/human/tracks/.{1,}", amazon.filenames)])
g_h.m.files <- c(g_h.m.files, amazon.filenames[grepl("data/human/methylation/.{1,}", amazon.filenames)])
i.temp <- amazon.filenames[grepl("data/human/hic/.{1,}", amazon.filenames)]
g_h.i.full <- c(g_h.i.full,  i.temp[!grepl(".rds", i.temp)])


# From 1-1,000,000-- ChIA-PET loops objects
if(length(g_h.c.full) != 0){
    g_h.c.names <- basename(file_path_sans_ext(g_h.c.full))
    g_h.c.list <- as.list(seq(1, length(g_h.c.names), by = 1) + 0)
    names(g_h.c.list) <- g_h.c.names
} else { g_h.c.list <- list(); g_h.c.full <- list()}

bigwig <- c(".bw", ".bigwig", ".bigWig")
# From 1,000,001-2,000,000-- ReadDepth Tracks-- bigwig
g_h.t.bw.full <- g_h.t.files[as.logical(rowSums(sapply(bigwig, grepl, g_h.t.files)))]
if(length(g_h.t.bw.full) != 0){
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
g_h.m.bw.full <- g_h.m.files[as.logical(rowSums(sapply(bigwig, grepl, g_h.m.files)))]
if(length(g_h.m.bw.full) != 0){
    g_h.m.bw.names <- basename(file_path_sans_ext(g_h.m.bw.full))
    g_h.m.bw.list <- as.list(seq(1, length(g_h.m.bw.names), by = 1) + 3000000)
    names(g_h.m.bw.list) <- g_h.m.bw.names
} else { g_h.m.bw.list <- list(); g_h.m.bw.full <- list()}

# From 4,000,001-5,000,000-- Methylation Tracks-- Bedgraph
g_h.m.bg.full <- g_h.m.files[grep(".bedgraph", g_h.m.files, fixed=T)]
if(length(g_h.m.bg.full) != 0){
    g_h.m.bg.names <- basename(file_path_sans_ext(g_h.m.bg.full))
    g_h.m.bg.list <- as.list(seq(1, length(g_h.m.bg.names), by = 1) + 4000000)
    names(g_h.m.bg.list) <- g_h.m.bg.names
} else { g_h.m.bg.list <- list(); g_h.m.bg.full <- list()}

# From 5,000,001-6,000,000-- HiC Tracks-- .rds
if(length(g_h.i.full) != 0){
    g_h.i.names <- basename(file_path_sans_ext(g_h.i.full))
    g_h.i.list <- as.list(seq(1, length(g_h.i.names), by = 1) + 5000000)
    names(g_h.i.list) <- g_h.i.names
} else {g_h.i.list <- list()}

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




## MOUSE INITIALIZATION ##
#Currently does not support HiC

# Import locally hosted data file names
g_m.c.full <- list.files("data/mouse/loops", full.names = TRUE)
g_m.t.files <- list.files("data/mouse/tracks", full.names = TRUE)
g_m.m.files <- list.files("data/mouse/methylation", full.names = TRUE)
g_m.i.full  <- list.files("data/mouse/hic", full.names = TRUE)

# Append amazon data
g_m.c.full <- c(g_m.c.full, amazon.filenames[grepl("data/mouse/loops/.{1,}", amazon.filenames)])
g_m.t.files <- c(g_m.t.files, amazon.filenames[grepl("data/mouse/tracks/.{1,}", amazon.filenames)])
g_m.m.files <- c(g_m.m.files, amazon.filenames[grepl("data/mouse/methylation/.{1,}", amazon.filenames)])
g_m.i.full <- c(g_m.i.full, amazon.filenames[grepl("data/mouse/hic/.{1,}", amazon.filenames)])


# From 1-1,000,000-- ChIA-PET loops objects
if(length(g_m.c.full) != 0){
    g_m.c.names <- basename(file_path_sans_ext(g_m.c.full))
    g_m.c.list <- as.list(seq(1, length(g_m.c.names), by = 1) + 0)
    names(g_m.c.list) <- g_m.c.names
} else { g_m.c.list <- list(); g_m.c.full <- list()}

bigwig <- c(".bw", ".bigwig", ".bigWig")
# From 1,000,001-2,000,000-- ReadDepth Tracks-- bigwig
g_m.t.bw.full <- g_m.t.files[as.logical(rowSums(sapply(bigwig, grepl, g_m.t.files)))]
if(length(g_m.t.bw.full) != 0){
    g_m.t.bw.names <- basename(file_path_sans_ext(g_m.t.bw.full))
    g_m.t.bw.list <- as.list(seq(1, length(g_m.t.bw.names), by = 1) + 1000000)
    names(g_m.t.bw.list) <- g_m.t.bw.names
} else { g_m.t.bw.list <- list(); g_m.t.bw.full <- list()}

# From 2,000,001-3,000,000-- ReadDepth Tracks-- Bedgraph
g_m.t.bg.full <- g_m.t.files[grep(".bedgraph", g_m.t.files, fixed=T)]
if(length(g_m.t.bg.full) != 0){
    g_m.t.bg.names <- basename(file_path_sans_ext(g_m.t.bg.full))
    g_m.t.bg.list <- as.list(seq(1, length(g_m.t.bg.names), by = 1) + 2000000)
    names(g_m.t.bg.list) <- g_m.t.bg.names
} else { g_m.t.bg.list <- list(); g_m.t.bg.full <- list() }

# From 3,000,001-4,000,000-- Methylation Tracks-- bigwig
# workaround for no methylation bigwigs
g_m.m.bw.full <- list()
temp <- sapply(bigwig, grepl, g_m.m.files)
if(!is.null(dim(temp))) g_m.m.bw.full <- try(g_m.m.files[as.logical(rowSums(temp))])
if(length(g_m.m.bw.full) != 0){
    g_m.m.bw.names <- basename(file_path_sans_ext(g_m.m.bw.full))
    g_m.m.bw.list <- as.list(seq(1, length(g_m.m.bw.names), by = 1) + 3000000)
    names(g_m.m.bw.list) <- g_m.m.bw.names
} else { g_m.m.bw.list <- list(); g_m.m.bw.full <- list()}

# From 4,000,001-5,000,000-- Methylation Tracks-- Bedgraph
g_m.m.bg.full <- g_m.m.files[grep(".bedgraph", g_m.m.files, fixed=T)]
if(length(g_m.m.bg.full) != 0){
    g_m.m.bg.names <- basename(file_path_sans_ext(g_m.m.bg.full))
    g_m.m.bg.list <- as.list(seq(1, length(g_m.m.bg.names), by = 1) + 4000000)
    names(g_m.m.bg.list) <- g_m.m.bg.names
} else { g_m.m.bg.list <- list(); g_m.m.bg.full <- list()}

# From 5,000,001-6,000,000-- HiC Tracks-- .rds
if(length(g_m.i.full) != 0){
    g_m.i.names <- basename(file_path_sans_ext(g_m.i.full))
    g_m.i.list <- as.list(seq(1, length(g_m.i.names), by = 1) + 5000000)
    names(g_m.i.list) <- g_h.i.names
} else {g_m.i.list <- list()}


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

