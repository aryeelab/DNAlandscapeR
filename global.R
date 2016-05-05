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

# From 0-1,000-- ChIA-PET loops objects
c.full <- list.files("data/loops", full.names = TRUE)
c.names <- basename(file_path_sans_ext(c.full))
c.list <- as.list(seq(1, length(c.names), by = 1) + 0)
names(c.list) <- c.names

# From 1,000-2,000-- Epigenetics Tracks
e.full <- list.files("data/tracks", full.names = TRUE)
e.names <- basename(file_path_sans_ext(e.full))
e.list <- as.list(seq(1, length(e.names), by = 1) + 1000)
names(e.list) <- e.names

# From 2,000-3,000-- DNA Methylation
m.full <- list.files("data/methylation", full.names = TRUE)
m.names <- basename(file_path_sans_ext(m.full))
m.list <- as.list(seq(1, length(m.names), by = 1) + 2000)
names(m.list) <- m.names

f.list <- append(c.list, append(e.list, m.list))
