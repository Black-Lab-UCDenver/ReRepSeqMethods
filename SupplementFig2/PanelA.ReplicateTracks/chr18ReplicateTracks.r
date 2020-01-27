library(GenomicRanges)
library(Gviz)
library(rtracklayer)
library(readr)
library(dplyr)



setwd("/path/to/files")

zero <- rtracklayer::import.bedGraph("noc3t0.chromsWeWantAndX.sort.bl.zerosFilled.1kbBin.1mbSmooth.bedgraph")
ten <- rtracklayer::import.bedGraph("noc3t10.chromsWeWantAndX.sort.bl.zerosFilled.1kbBin.1mbSmooth.bedgraph")
fifteen <- rtracklayer::import.bedGraph("noc3t15.chromsWeWantAndX.sort.bl.zerosFilled.1kbBin.1mbSmooth.bedgraph")
twentyFive <- rtracklayer::import.bedGraph("noc3t25.chromsWeWantAndX.sort.bl.zerosFilled.1kbBin.1mbSmooth.bedgraph")


ERD <- rtracklayer::import.bed("merge350kERD.bed")
LRD <- rtracklayer::import.bed("merge350kLRD.bed") 

erdHighLight <- read_tsv("merge350kERD.bed", col_names = FALSE)


chromSize <- read_tsv("hg38.chrom.sizes.sort", col_names = FALSE)

  
i <- "chr18"

  thisChromHighlights <- filter(erdHighLight, X1==i)
  thisChromHighlights$X4 <- thisChromHighlights$X3 - thisChromHighlights$X2
  
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
  
  

  
  thisChromSize <- filter(chromSize, X1==paste(i))
  tickSpacing <- thisChromSize$X2/4
  first <- as.numeric((round(tickSpacing/1000000))*1000000)
  second <- as.numeric((round((tickSpacing*2)/1000000))*1000000)
  third <- as.numeric((round((tickSpacing*3)/1000000))*1000000)
  fourth <- as.numeric((round((tickSpacing*4)/1000000))*1000000)
  fifth <- as.numeric((round((tickSpacing*5)/1000000))*1000000)
  sixth <- as.numeric((round((tickSpacing*6)/1000000))*1000000)
  
  gtrack <- GenomeAxisTrack(fontcolor = "black", col="black", labelPos="below")
  
  highlight <- HighlightTrack(trackList = list(zeroTrack, tenTrack, fifteenTrack, twentyFiveTrack), start= c(thisChromHighlights$X2), width=c(thisChromHighlights$X4), chromosome = i, col="#9ebeff", fill="#9ebeff", inBackground=TRUE)
  
  tiff(paste0("rep3.", chromName, ".tiff"), width = 325, height = 425)
  outPlot<- plotTracks(list(highlight, bedERD, bedLRD, gtrack), type = "hist", fontsize=12, ticksAt=c(first, second, third), cex=1, cex.axis=1, fontfamily="sans")
  dev.off()
  


