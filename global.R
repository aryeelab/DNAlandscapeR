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
library(shinyFiles)
library(DT)
options(warn=-1)

# Import locally hosted data file names
c.full <- list.files("data/loops", full.names = TRUE)
t.files <- list.files("data/tracks", full.names = TRUE)
m.files <- list.files("data/methylation/", full.names = TRUE)

# From 1-1,000-- ChIA-PET loops objects
if(length(c.full) != 0){
    c.names <- basename(file_path_sans_ext(c.full))
    c.list <- as.list(seq(1, length(c.names), by = 1) + 0)
    names(c.list) <- c.names
} else { c.list <- list(); c.full <- list()}

# From 1,001-2,000-- ReadDepth Tracks-- bigwig
t.bw.full <- t.files[grep(".bw", t.files, fixed=T)]
if(length(t.bw.full) != 0){
    t.bw.names <- basename(file_path_sans_ext(t.bw.full))
    t.bw.list <- as.list(seq(1, length(t.bw.names), by = 1) + 1000)
    names(t.bw.list) <- t.bw.names
} else { t.bw.list <- list(); t.bw.full <- list()}

# From 2,001-3,000-- ReadDepth Tracks-- Bedgraph
t.bg.full <- t.files[grep(".bedgraph", t.files, fixed=T)]
if(length(t.bg.full) != 0){
    t.bg.names <- basename(file_path_sans_ext(t.bg.full))
    t.bg.list <- as.list(seq(1, length(t.bg.names), by = 1) + 2000)
    names(t.bg.list) <- t.bg.names
} else { t.bg.list <- list(); t.bg.full <- list() }

# From 3,001-4,000-- Methylation Tracks-- bigwig
m.bw.full <- m.files[grep(".bw", m.files, fixed=T)]
if(length(m.bw.full) != 0){
    m.bw.names <- basename(file_path_sans_ext(m.bw.full))
    m.bw.list <- as.list(seq(1, length(m.bw.names), by = 1) + 3000)
    names(m.bw.list) <- m.bw.names
} else { m.bw.list <- list(); m.bw.full <- list()}

# From 4,001-5,000-- Methylation Tracks-- Bedgraph
m.bg.full <- m.files[grep(".bedgraph", m.files, fixed=T)]
if(length(m.bg.full) != 0){
    m.bg.names <- basename(file_path_sans_ext(m.bg.full))
    m.bg.list <- as.list(seq(1, length(m.bg.names), by = 1) + 4000)
    names(m.bg.list) <- m.bg.names
} else { m.bg.list <- list(); m.bg.full <- list()}

f.list <- append(c.list, append(append(t.bw.list, t.bg.list), append(m.bw.list, m.bg.list)))
