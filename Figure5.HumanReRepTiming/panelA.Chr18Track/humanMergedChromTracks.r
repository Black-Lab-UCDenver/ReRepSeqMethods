######################ReadMe####################
#This code will produce the genomic tracks presented in figure 5 panel a. The input data is supplied in this directory.
# When attempting to run this script, make sure to set the working directory to location in which the input files are saved.
# The input files are consist of merged, normalized and scaled, bedgraphs from the human noc release experiment. Specifally, 
# time points 0, 10 hours, 15 hours, 25 hours, and two cell cycle samples. 
#
# Additionally, ERD and LRD files are suplied for annotations. These ERD and LRD
# files are a merged form of the original files human ERD and LRD files. When attempting top make smaller figures, Gviz will
# generate additional lines of annotation if the white space between two freatures of the same annotation is not small enough to display.
# This creates inconcsistently sized figures that adds little information to the graphic. We generated merged ERD or LRD regions into
# into a single regions that were seperated by 350kb or less. This allowed us to display all ERDs and LRDs over a whole chromosome
# without distoring the figure. 
#
#The ERD regions are loaded a second time as a highlighted track to help visualize thier 
# boundaries with respect to signal in the tracks.
#
#This script will produce a chromosome plot for every chromsome listed in the bedgraph files.
############################################


#Load libraries
library(GenomicRanges)
library(Gviz)
library(rtracklayer)
library(readr)
library(dplyr)


#Change this directory to location of the downloaded inout files.
setwd("/path/to/files")

#Load the input files for each timepoint and the 2 cell cycle control.
  zero <- rtracklayer::import.bedGraph("noct0.mergedAve.bedgraph")
  ten <- rtracklayer::import.bedGraph("noct10.mergedAve.bedgraph")
  fifteen <- rtracklayer::import.bedGraph("noct15.mergedAve.bedgraph")
  twentyFive <- rtracklayer::import.bedGraph("noct25.mergedAve.bedgraph")
  twoCellCycle <- rtracklayer::import.bedGraph("h2cc.mergedAve.bedgraph")

#load erd and lrd files  
  ERD <- rtracklayer::import.bed("ERD.bed")
  LRD <- rtracklayer::import.bed("LRD.bed") 

#load the red regionns as a highlight track  
  erdHighLight <- read_tsv("ERD.bed", col_names = FALSE)

#Make a list of chromosomes based on the chromosomes found in the input bedgraphs    
  chrList <- unique(zero@seqinfo@seqnames)

  chromSize <- read_tsv("hg38.chrom.sizes.sort", col_names = FALSE)
  for(i in chrList){
    
    #Define the boudaries for each ERD region to highlight in the tracks.
    thisChromHighlights <- filter(erdHighLight, X1==i)
    thisChromHighlights$X4 <- thisChromHighlights$X3 - thisChromHighlights$X2
  
    
    # generate a track for each time point and annotation file.  
    chromName <- i
    zeroTrack <- DataTrack(zero, genome = "hg38", 
                        type = "hist", 
                        chromosome = i, 
                        name = "0hr", 
                        background.title = "white", 
                        col.title = "black", 
                        col.axis = "black", 
                        col.histogram = "blue", 
                        ylim = c(0,300))
    tenTrack <- DataTrack(ten, genome = "hg38", 
                      type = "hist", 
                      chromosome = i, 
                      name = "10hr", 
                      background.title = "white", 
                      col.title = "black", 
                      col.axis = "black", 
                      col.histogram = "blue", 
                      ylim = c(0,300))
    fifteenTrack <- DataTrack(fifteen, genome = "hg38", 
                      type = "hist", 
                      chromosome = i, 
                      name = "15hr", 
                      background.title = "white", 
                      col.title = "black", 
                      col.axis = "black", 
                      col.histogram = "blue", 
                      ylim = c(0,300))
    twentyFiveTrack <- DataTrack(twentyFive, genome = "hg38", 
                      type = "hist", 
                      chromosome = i, 
                      name = "25hr", 
                      background.title = "white", 
                      col.title = "black", 
                      col.axis = "black", 
                      col.histogram = "blue", 
                      ylim = c(0,300))
    
    twoCellCycleTrack <- DataTrack(twoCellCycle, genome = "hg38", 
                                 type = "hist", 
                                 chromosome = i, 
                                 name = "2CC", 
                                 background.title = "white", 
                                 col.title = "black", 
                                 col.axis = "black", 
                                 col.histogram = "black", 
                                 ylim = c(0,300))
    
    bedERD <- AnnotationTrack(ERD, genome = "hg38", 
                              chromosome = i,
                              # name = "ERD",
                              col = "blue",
                              fill = "blue",
                              collapse = TRUE, 
                              background.panel="white", 
                              background.title="white", 
                              stacking="dense", 
                              min.height=1)
    
    bedLRD <- AnnotationTrack(LRD, genome = "hg38", 
                              chromosome = i,
                              # name = "LRD",
                              col = "red",
                              fill = "red",
                              collapse = TRUE, 
                              background.panel="white", 
                              background.title="white", 
                              stacking="dense", 
                              min.height=1)
    
  
    #Dymaically position chromosome tick marks proportionally to the size of the chromosome
    thisChromSize <- filter(chromSize, X1==paste(i))
    tickSpacing <- thisChromSize$X2/4
    first <- as.numeric((round(tickSpacing/1000000))*1000000)
    second <- as.numeric((round((tickSpacing*2)/1000000))*1000000)
    third <- as.numeric((round((tickSpacing*3)/1000000))*1000000)
    fourth <- as.numeric((round((tickSpacing*4)/1000000))*1000000)
    fifth <- as.numeric((round((tickSpacing*5)/1000000))*1000000)
    sixth <- as.numeric((round((tickSpacing*6)/1000000))*1000000)
    
    #add a chromsome axis
    gtrack <- GenomeAxisTrack(fontcolor = "black", col="black", labelPos="below")
    
    # Choose which tracks to highlight
    highlight <- HighlightTrack(trackList = list(zeroTrack, tenTrack, fifteenTrack, twentyFiveTrack, twoCellCycleTrack), start= c(thisChromHighlights$X2), width=c(thisChromHighlights$X4), chromosome = i, col="#9ebeff", fill="#9ebeff", inBackground=TRUE)
    
    #save the outputs.
    tiff(paste0("merged.", chromName, ".tiff"), width = 325, height = 425)
    outPlot<- plotTracks(list(highlight, bedERD, bedLRD, gtrack), type = "hist", fontsize=12, ticksAt=c(first, second, third), cex=1, cex.axis=1, fontfamily="sans")
    dev.off()
    
  }
  
