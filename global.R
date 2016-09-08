# Global variables for both server.R and ui.R to reference

source("www/libraries.R")
source("www/adv-shiny.R")
source("www/addbucket.R")
source("www/js.R")
source("www/libraryNormHiC.R")

options(shiny.error=browser)
options(shiny.maxRequestSize=1*1024^3) #1 GB Max file size
options(warn=-1)
options(scipen=999)

# Get commit ID
sha <- git2r::commits()[[1]]@sha
short_sha <- substr(sha, 1, 7)

# Set up the local S3 cache dir
cache_dir <- tempdir()

default_chr <- "9"
default_start <- 21912689
default_end <- 22216233
ucsc_coord <- paste0("chr", default_chr, ":", as.character(default_start), "-", as.character(default_end))

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
t <- unlist(get_bucket(bucket = "dnalandscaper"))
amazon.filenames <- paste(amazon, t[grep("data", t)], sep = "/")

## HUMAN INITIALIZATION ##

# Import locally hosted data file names
g_h.c.full <- list.files("data/human/loops", full.names = TRUE)
g_h.t.files <- list.files("data/human/tracks", full.names = TRUE)
g_h.m.files <- list.files("data/human/methylation", full.names = TRUE)

# Local HiC-- voided currently
g_h.i.samples <- list.files("data/human/hic", full.names = FALSE, recursive = FALSE)
res.temp <- list.dirs(paste0("data/human/hic/", g_h.i.samples), full.names = FALSE, recursive = FALSE)
g_h.i.res <- lapply(g_h.i.samples, function(t){unlist(strsplit(res.temp[grepl(t, res.temp)],split="_"))[c(FALSE,TRUE)]})
g_h.i.full <- list.files("data/human/hic", recursive = TRUE, full.names = TRUE)

# Append Amazon data
g_h.c.full <- c(g_h.t.files, amazon.filenames[grepl("data/human/loops/.{1,}", amazon.filenames)])
g_h.t.files <- c(g_h.t.files, amazon.filenames[grepl("data/human/tracks/.{1,}", amazon.filenames)])
g_h.m.files <- c(g_h.m.files, amazon.filenames[grepl("data/human/methylation/.{1,}", amazon.filenames)])

# Append Amazon HiC data
i.temp <- amazon.filenames[grepl("data/human/hic/.{1,}", amazon.filenames)]
i.base <- basename(i.temp)
amazon.hic.samples <- file_path_sans_ext(file_path_sans_ext(i.base[grepl(".sparseHiC.meta", i.base)]))
g_h.i.samples <- c(g_h.i.samples, amazon.hic.samples)
g_h.i.res <- c(g_h.i.res, lapply(amazon.hic.samples, function(t){
    restab <- read.table(i.temp[grepl(t, i.temp) & grepl(".sparseHiC.meta", i.temp)], skip = 3)
    strsplit(as.character(restab[1,1]), split = ",")[[1]]
}))
names(g_h.i.res) <- g_h.i.samples
g_h.i.full <- c(g_h.i.full, i.temp[grepl(".rds", i.temp)])


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
if(length(g_h.i.samples) != 0){
    g_h.i.list <- as.list(seq(1, length(g_h.i.samples), by = 1) + 5000000)
    names(g_h.i.list) <- g_h.i.samples
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

# Import locally hosted data file names
g_m.c.full <- list.files("data/mouse/loops", full.names = TRUE)
g_m.t.files <- list.files("data/mouse/tracks", full.names = TRUE)
g_m.m.files <- list.files("data/mouse/methylation", full.names = TRUE)

# Local HiC data
g_m.i.samples <- list.dirs("data/mouse/hic/", full.names = FALSE, recursive = FALSE)
m.res.temp <- list.dirs(paste0("data/mouse/hic/", g_m.i.samples), full.names = FALSE, recursive = FALSE)
g_m.i.res <- lapply(g_m.i.samples, function(t){unlist(strsplit(m.res.temp[grepl(t, m.res.temp)],split="_"))[c(FALSE,TRUE)]})
g_m.i.full <- list.files("data/mouse/hic", recursive = TRUE, full.names = TRUE)

# Append amazon data
g_m.c.full <- c(g_m.c.full, amazon.filenames[grepl("data/mouse/loops/.{1,}", amazon.filenames)])
g_m.t.files <- c(g_m.t.files, amazon.filenames[grepl("data/mouse/tracks/.{1,}", amazon.filenames)])
g_m.m.files <- c(g_m.m.files, amazon.filenames[grepl("data/mouse/methylation/.{1,}", amazon.filenames)])

# Append Amazon HiC data
i.temp <- amazon.filenames[grepl("data/mouse/hic/.{1,}", amazon.filenames)]
i.base <- basename(i.temp)
amazon.hic.samples <- file_path_sans_ext(file_path_sans_ext(i.base[grepl(".sparseHiC.meta", i.base)]))
g_m.i.samples <- c(g_m.i.samples, amazon.hic.samples)
g_m.i.res <- c(g_m.i.res, lapply(amazon.hic.samples, function(t){
    restab <- read.table(i.temp[grepl(t, i.temp) & grepl(".sparseHiC.meta", i.temp)], skip = 3)
    strsplit(as.character(restab[1,1]), split = ",")[[1]]
}))
names(g_m.i.res) <- g_m.i.samples
g_m.i.full <- c(g_m.i.full, i.temp[grepl(".rds", i.temp)])

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
temp <- matrix(sapply(bigwig, grepl, g_m.m.files), ncol = 3)
g_m.m.bw.full <- g_m.m.files[as.logical(rowSums(temp))]
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
if(length(g_m.i.samples) != 0){
    g_m.i.list <- as.list(seq(1, length(g_m.i.samples), by = 1) + 5000000)
    names(g_m.i.list) <- g_m.i.samples
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

