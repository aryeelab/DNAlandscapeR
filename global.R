# Global variables for both server.R and ui.R to reference
library(tools)
library(shiny)
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

uploadchoices <- list("Loops/.rds" = 1,
                      "ReadDepth/.bigWig" = 2,
                      "ReadDepth/.bedgraph" = 3,
                      "Methylation/.bigWig" = 4,
                      "Methylation/.bedgraph" = 5)

# Import locally hosted data file names
c.full <- list.files("data/loops", full.names = TRUE)
t.files <- list.files("data/tracks", full.names = TRUE)
m.files <- list.files("data/methylation", full.names = TRUE)

# From 1-1,000,000-- ChIA-PET loops objects
if(length(c.full) != 0){
    c.names <- basename(file_path_sans_ext(c.full))
    c.list <- as.list(seq(1, length(c.names), by = 1) + 0)
    names(c.list) <- c.names
} else { c.list <- list(); c.full <- list()}

bigwig <- c(".bw", ".bigwig", ".bigWig")
# From 1,000,001-2,000,000-- ReadDepth Tracks-- bigwig
t.bw.full <- t.files[as.logical(rowSums(sapply(bigwig, grepl, t.files)))]
if(length(t.bw.full) != 0){
    t.bw.names <- basename(file_path_sans_ext(t.bw.full))
    t.bw.list <- as.list(seq(1, length(t.bw.names), by = 1) + 1000000)
    names(t.bw.list) <- t.bw.names
} else { t.bw.list <- list(); t.bw.full <- list()}

# From 2,000,001-3,000,000-- ReadDepth Tracks-- Bedgraph
t.bg.full <- t.files[grep(".bedgraph", t.files, fixed=T)]
if(length(t.bg.full) != 0){
    t.bg.names <- basename(file_path_sans_ext(t.bg.full))
    t.bg.list <- as.list(seq(1, length(t.bg.names), by = 1) + 2000000)
    names(t.bg.list) <- t.bg.names
} else { t.bg.list <- list(); t.bg.full <- list() }

# From 3,000,001-4,000,000-- Methylation Tracks-- bigwig
m.bw.full <- m.files[as.logical(rowSums(sapply(bigwig, grepl, m.files)))]
if(length(m.bw.full) != 0){
    m.bw.names <- basename(file_path_sans_ext(m.bw.full))
    m.bw.list <- as.list(seq(1, length(m.bw.names), by = 1) + 3000000)
    names(m.bw.list) <- m.bw.names
} else { m.bw.list <- list(); m.bw.full <- list()}

# From 4,000,001-5,000,000-- Methylation Tracks-- Bedgraph
m.bg.full <- m.files[grep(".bedgraph", m.files, fixed=T)]
if(length(m.bg.full) != 0){
    m.bg.names <- basename(file_path_sans_ext(m.bg.full))
    m.bg.list <- as.list(seq(1, length(m.bg.names), by = 1) + 4000000)
    names(m.bg.list) <- m.bg.names
} else { m.bg.list <- list(); m.bg.full <- list()}

f.list <- append(c.list, append(append(t.bw.list, t.bg.list), append(m.bw.list, m.bg.list)))

# Now order it
d <- data.frame(unlist(f.list))
d$names <- rownames(d)
d <- d[order(rownames(d)), ]
f.list <- list()
for(k in 1:dim(d)[1]){
    x <- list(d[k,1])
    names(x) <- (d[k,2])
    f.list <- append(f.list, x)
}