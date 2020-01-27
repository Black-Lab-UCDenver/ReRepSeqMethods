###########################ReadMe#############
#This code was used to generate the chromosome tracks featured in figure 3 panel a. 
#Notably, as written, this will output figure panels for all chromosomes.
#All supporting input files are privded on github: https://github.com/blacklabUCD/ReRepSeqMethods
#This code will take as input the following files:
#1) ERD highlight tracks
#2) ERD regions, which is the same file as the ERD highlight tracks, but is imported twice to fulfill two different functions.
#3) LRD regions
#4) Confirmed ARS
#5) Topological domains
#6) merged bedgraphs for the yeast replication timing experiment
#7) chromosome sizes.
# Please note this section of code relies on the instalation of R3.6 to run Gviz.


#Load libraries
library(GenomicRanges)
library(Gviz)
library(rtracklayer)
library(readr)
library(dplyr)

#set the working directory to the location of the files
setwd("/path/to/files")

#Import erd/lrd/tad/ars bed files
erdHighLight <- read_tsv("erdYeast.merged.10k.bed", col_names = FALSE)
ERD <- rtracklayer::import.bed("erdYeast.merged.10k.bed")
LRD <- rtracklayer::import.bed("lrdYeast.merged.10k.bed")
ARS <- rtracklayer::import.bed("arsConfirmed.bed")
TAD <- rtracklayer::import.bed("tad.bed")


#Import samples to be analyzed
inSample0 <- rtracklayer::import.bedGraph("0min.merged.final.ave.bedgraph")
inSample15 <- rtracklayer::import.bedGraph("15min.merged.final.ave.bedgraph")
inSample20 <- rtracklayer::import.bedGraph("20min.merged.final.ave.bedgraph")
inSample25 <- rtracklayer::import.bedGraph("25min.merged.final.ave.bedgraph")
inSample30 <- rtracklayer::import.bedGraph("30min.merged.final.ave.bedgraph")
inSample35 <- rtracklayer::import.bedGraph("35min.merged.final.ave.bedgraph")
inSample60 <- rtracklayer::import.bedGraph("60min.merged.final.ave.bedgraph")
inSample80 <- rtracklayer::import.bedGraph("80min.merged.final.ave.bedgraph")
inSampley3cc <- rtracklayer::import.bedGraph("y3cc.merged.final.ave.bedgraph")


#Make a list of all chromosomes found in the samples.
chrList <- unique(inSample0@seqinfo@seqnames)

#Import chromsome sizes for reference.
chromSize <- read_tsv("sacCer3.chrom.sizes", col_names = FALSE)

#This is the function that loops through all the chromosomes in the chrList file to generate all the tracks.
for(i in chrList){
  
  # "i" is defined as an individual chromosome. Using this variable, all necessary information for the chromosome
  # corresponding to "i" is retreuved.
  
  #get "i" erd highlights
  thisChromHighlights <- filter(erdHighLight, X1==i)
  
  #calculate the width of each erd region to be highlighted
  thisChromHighlights$X4 <- thisChromHighlights$X3 - thisChromHighlights$X2
  
  #Get the name of the chromosome, which will be used later when saving the file
  chromName <- i
  
  #Process the data for each timepoint and annotation to be plotted
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
  
  y3ccTrack <- DataTrack(inSampley3cc, genome = "sacCer3", 
                     type = "hist", 
                     chromosome = i, 
                     name = i, 
                     background.title = "white", 
                     col.title = "black", 
                     col.axis = "black", 
                     col.histogram = "black", 
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
  
  
 
 #Using the chromosome size file, divide the current chromosome into equal sections to generate evenly spaced chromosome arm labels.
  thisChromSize <- filter(chromSize, X1==paste(i))
  tickSpacing <- thisChromSize$X2/4
  first <- as.numeric((round(tickSpacing/1000))*1000)
  second <- as.numeric((round((tickSpacing*2)/1000))*1000)
  third <- as.numeric((round((tickSpacing*3)/1000))*1000)
  fourth <- as.numeric((round((tickSpacing*4)/1000))*1000)
  fifth <- as.numeric((round((tickSpacing*5)/1000))*1000)
  sixth <- as.numeric((round((tickSpacing*6)/1000))*1000)
  
  # Add a genomic axis to the plot
  gtrack <- GenomeAxisTrack(fontcolor = "black", col="black", labelPos="below")
  
 # Determine which tracks arfe to be hihglighted and in which order.
  highlight <- HighlightTrack(trackList = list(min0, min15, min20, min25, min30, min35, min60, min80, y3ccTrack), start= c(thisChromHighlights$X2), width=c(thisChromHighlights$X4), chromosome = i, col="#9ebeff", fill="#9ebeff", inBackground=TRUE)
  
  # Automatically save the resulting graphic to the current directory and make it after the chomosome it was derived.
  tiff(paste0("./", chromName, ".tiff"), width = 300, height = 400)
  outPlot<- plotTracks(list(highlight, bedARS, bedERD, bedLRD, bedTAD, gtrack), type = "hist", fontsize=12, ticksAt=c(first, second, third), cex=1, cex.axis=1, fontfamily="sans")
  dev.off()
  
}

