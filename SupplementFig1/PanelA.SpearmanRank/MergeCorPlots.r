library(readr)
library(dplyr)
library(corrplot)
library(RColorBrewer)

setwd("/path/to/data")

inFilePaths <- list.files(path = ".", pattern = glob2rx("*min.bed"), full.names = TRUE)

for(inFilePath in inFilePaths){
  
  rm(dfIn, dfNext, name, dfSplit, aveOut, aveFinal)
  
  dfIn <- read_tsv(inFilePath, col_names = FALSE)
  name <- gsub(".bed", "", inFilePath)
  name <- gsub("./", "", name)
  
  dfSplit<- split(dfIn, dfIn$X8)
  
  for(i in dfSplit){
    ars <- i[1,8]
    numRegions <- length(i$X1)
    aveVal <- sum(i$X4)/numRegions
   
    outAve <- as.data.frame(c(ars, aveVal))
    colnames(outAve) <- c("Gene", name) 
    
    if(exists("aveFinal") == FALSE){
      aveFinal <- outAve
    }else{
      aveFinal <- rbind(aveFinal, outAve)
    }
  }
  

  if(exists("dfOut")==FALSE){
    dfOut <- aveFinal

  }else{
    dfOut<- full_join(dfOut, aveFinal, by=c("Gene"="Gene"))
  }
  
  
}

dfOut[is.na(dfOut)] <- 0


dfCor <- select(dfOut, s1.0min, s2.0min, s3.0min, s1.15min, s2.15min, s3.15min, s1.20min, s2.20min, s3.20min, s1.25min, s2.25min, s3.25min, s1.30min, s2.30min, s3.30min, s1.35min, s2.35min, s3.35min, s1.60min, s2.60min, s3.60min, s1.80min, s2.80min, s3.80min)

m<- cor(dfCor, method="spearman")
dfM <- as.data.frame(m)
write_tsv(dfM, "CorMatrixYeastRepTimingSpearman.txt", col_names = TRUE)



