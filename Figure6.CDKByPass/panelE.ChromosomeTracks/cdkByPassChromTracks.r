
#load libraries
library(GenomicRanges)
library(Gviz)
library(rtracklayer)


# set working directory
setwd("/path/to/files")

#load input files
  ERD <- rtracklayer::import.bed("/media/phil/external2/erdYeast.merged.10k.bed")
  LRD <- rtracklayer::import.bed("/media/phil/external2/lrdYeast.merged.10k.bed")
  ARS <- rtracklayer::import.bed("arsConfirmed.bed")
  TAD <- rtracklayer::import.bed("tad.bed")
  
  inSampleGal <- rtracklayer::import.bedGraph("/media/phil/external2/yeastFinal/cdkbypass/g.mergedAve.100bin.10kbsmooth.bedgraph")
  inSampleRaf <- rtracklayer::import.bedGraph("/media/phil/external2/yeastFinal/cdkbypass/r.mergedAve.100bin.10kbsmooth.bedgraph")
  chrList <- unique(inSampleGal@seqinfo@seqnames)

  chromSize <- read_tsv("sacCer3.chrom.sizes", col_names = FALSE)

#make an annotated track for each chromosome   
    for(i in chrList){
    
    chromName <- i
    Gal <- DataTrack(inSampleGal, genome = "sacCer3", 
                        type = "hist", 
                        chromosome = i, 
                        name = i, 
                        background.title = "white", 
                        col.title = "black", 
                        col.axis = "black", 
                        col.histogram = "red", 
                        ylim = c(0,10))
    Raf <- DataTrack(inSampleRaf, genome = "sacCer3", 
                     type = "hist", 
                     chromosome = i, 
                     name = i, 
                     background.title = "white", 
                     col.title = "black", 
                     col.axis = "black", 
                     col.histogram = "black", 
                     ylim = c(0,10))
    
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
    
    
  # dynamically define the location of the chromosome labels based on chromosome size
    thisChromSize <- filter(chromSize, X1==paste(i))
    tickSpacing <- thisChromSize$X2/4
    first <- as.numeric((round(tickSpacing/1000))*1000)
    second <- as.numeric((round((tickSpacing*2)/1000))*1000)
    third <- as.numeric((round((tickSpacing*3)/1000))*1000)
    fourth <- as.numeric((round((tickSpacing*4)/1000))*1000)
    fifth <- as.numeric((round((tickSpacing*5)/1000))*1000)
    sixth <- as.numeric((round((tickSpacing*6)/1000))*1000)
    
    gtrack <- GenomeAxisTrack(fontcolor = "black", col="black", labelPos="below")
    
    #save the files
    tiff(paste0("cdkByPass.Gal.Raf.", chromName, ".tiff"), width = 300, height = 200)
    outPlot<- plotTracks(list(Gal, Raf, bedARS, bedERD, bedLRD, bedTAD, gtrack), type = "hist", fontsize=12, ticksAt=c(first, second, third), cex=1, cex.axis=1, fontfamily="sans")
    dev.off()
    
  }
  
