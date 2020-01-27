
#load libraries
library(GenomicRanges)
library(Gviz)
library(rtracklayer)
library(readr)
library(dplyr)


#set the working directory
setwd("/path/to/files")

erdHighLight <- read_tsv("erdYeast.merged.10k.bed", col_names = FALSE)

ERD <- rtracklayer::import.bed("erdYeast.merged.10k.bed")
LRD <- rtracklayer::import.bed("lrdYeast.merged.10k.bed")
ARS <- rtracklayer::import.bed("arsConfirmed.bed")
TAD <- rtracklayer::import.bed("tad.bed")

inSample0 <- rtracklayer::import.bedGraph("s3.0min.sort.blRpmMito100bin10kbSmooth.bg")
inSample15 <- rtracklayer::import.bedGraph("s3.15min.sort.blRpmMito100bin10kbSmooth.bg")
inSample20 <- rtracklayer::import.bedGraph("s3.20min.sort.blRpmMito100bin10kbSmooth.bg")
inSample25 <- rtracklayer::import.bedGraph("s3.25min.sort.blRpmMito100bin10kbSmooth.bg")
inSample30 <- rtracklayer::import.bedGraph("s3.30min.sort.blRpmMito100bin10kbSmooth.bg")
inSample35 <- rtracklayer::import.bedGraph("s3.35min.sort.blRpmMito100bin10kbSmooth.bg")
inSample60 <- rtracklayer::import.bedGraph("s3.60min.sort.blRpmMito100bin10kbSmooth.bg")
inSample80 <- rtracklayer::import.bedGraph("s3.80min.sort.blRpmMito100bin10kbSmooth.bg")

chromSize <- read_tsv("sacCer3.chrom.sizes", col_names = FALSE)

i <- "chrII"  
  thisChromHighlights <- filter(erdHighLight, X1==i)
  thisChromHighlights$X4 <- thisChromHighlights$X3 - thisChromHighlights$X2
  chromName <- i
  min0 <- DataTrack(inSample0, genome = "sacCer3", 
                    type = "hist", 
                    chromosome = i, 
                    name = i, 
                    background.title = "white", 
                    col.title = "black", 
                    col.axis = "black", 
                    col.histogram = "blue", 
                    ylim = c(0,30))
  
  min15 <- DataTrack(inSample15, genome = "sacCer3", 
                     type = "hist", 
                     chromosome = i, 
                     name = i, 
                     background.title = "white", 
                     col.title = "black", 
                     col.axis = "black", 
                     col.histogram = "blue", 
                     ylim = c(0,30))
  
  min20 <- DataTrack(inSample20, genome = "sacCer3", 
                     type = "hist", 
                     chromosome = i, 
                     name = i, 
                     background.title = "white", 
                     col.title = "black", 
                     col.axis = "black", 
                     col.histogram = "blue", 
                     ylim = c(0,30))
  
  min25 <- DataTrack(inSample25, genome = "sacCer3", 
                     type = "hist", 
                     chromosome = i, 
                     name = i, 
                     background.title = "white", 
                     col.title = "black", 
                     col.axis = "black", 
                     col.histogram = "blue", 
                     ylim = c(0,30))
  
  min30 <- DataTrack(inSample30, genome = "sacCer3", 
                     type = "hist", 
                     chromosome = i, 
                     name = i, 
                     background.title = "white", 
                     col.title = "black", 
                     col.axis = "black", 
                     col.histogram = "blue", 
                     ylim = c(0,30))
  
  min35 <- DataTrack(inSample35, genome = "sacCer3", 
                     type = "hist", 
                     chromosome = i, 
                     name = i, 
                     background.title = "white", 
                     col.title = "black", 
                     col.axis = "black", 
                     col.histogram = "blue", 
                     ylim = c(0,30))
  
  min60 <- DataTrack(inSample60, genome = "sacCer3", 
                     type = "hist", 
                     chromosome = i, 
                     name = i, 
                     background.title = "white", 
                     col.title = "black", 
                     col.axis = "black", 
                     col.histogram = "blue", 
                     ylim = c(0,30))
  
  min80 <- DataTrack(inSample80, genome = "sacCer3", 
                     type = "hist", 
                     chromosome = i, 
                     name = i, 
                     background.title = "white", 
                     col.title = "black", 
                     col.axis = "black", 
                     col.histogram = "blue", 
                     ylim = c(0,30))
  
  
  bedERD <- AnnotationTrack(ERD, genome = "sacCer3", 
                            chromosome = i,
                            # name = "ERD",
                            col = "blue",
                            fill = "blue",
                            collapse = TRUE, 
                            background.panel="white", 
                            background.title="white")
  
  bedLRD <- AnnotationTrack(LRD, genome = "sacCer3", 
                            chromosome = i,
                            # name = "LRD",
                            col = "red",
                            fill = "red",
                            collapse = TRUE, 
                            background.panel="white", 
                            background.title="white")
  
  bedARS <- AnnotationTrack(ARS, genome = "sacCer3", 
                            chromosome = i,
                            # name = "LRD",
                            col = "darkgray",
                            fill = "darkgray", 
                            stacking = "dense", 
                            collapse = TRUE, 
                            background.panel="white", 
                            background.title="white")
  
  bedTAD <- AnnotationTrack(TAD, genome = "sacCer3", 
                            chromosome = i,
                            # name = "LRD",
                            col = "green",
                            fill = "green",
                            collapse = TRUE, 
                            background.panel="white", 
                            background.title="white")
  
  

  
  thisChromSize <- filter(chromSize, X1==paste(i))
  tickSpacing <- thisChromSize$X2/4
  first <- as.numeric((round(tickSpacing/1000))*1000)
  second <- as.numeric((round((tickSpacing*2)/1000))*1000)
  third <- as.numeric((round((tickSpacing*3)/1000))*1000)
  fourth <- as.numeric((round((tickSpacing*4)/1000))*1000)
  fifth <- as.numeric((round((tickSpacing*5)/1000))*1000)
  sixth <- as.numeric((round((tickSpacing*6)/1000))*1000)
  
  gtrack <- GenomeAxisTrack(fontcolor = "black", col="black", labelPos="below")
  
  
  highlight <- HighlightTrack(trackList = list(min0, min15, min20, min25, min30, min35, min60, min80), start= c(thisChromHighlights$X2), width=c(thisChromHighlights$X4), chromosome = i, col="#9ebeff", fill="#9ebeff", inBackground=TRUE)
  
  
  tiff(paste0("rep3.", chromName, ".tiff"), width = 300, height = 400)
  outPlot<- plotTracks(list(highlight, bedARS, bedERD, bedLRD, bedTAD, gtrack), type = "hist", fontsize=12, ticksAt=c(first, second, third), cex=1, cex.axis=1, fontfamily="sans")
  dev.off()
  

