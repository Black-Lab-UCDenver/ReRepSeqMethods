
#Load libraries
library(readr)
library(dplyr)
library(reshape2)
library(ggplot2)


#set working directory
setwd("/path/to/files")

#import left side of transition to timing domain midpoint
dfLeft <- read_tsv("noct15.100kbS.left.bed", col_names = FALSE)

#split by region
leftSplit <- split(dfLeft, dfLeft$X8)

#apply code to all regions in list  
for(i in leftSplit){ 
  dfw <- i
  lengthLeft <- round(length(dfw$X4)/100)
  
  #resize and scale regions to fit a standard window.  
  splitVecs<- split(dfw$X4, ceiling(seq_along(dfw$X4)/lengthLeft))
  
  for (n in splitVecs){
    meanN<- mean(n)
    if(exists("outAve")==TRUE){
      outAve<- c(outAve, meanN)
    }else{
      outAve <- meanN
    }
  }
  
  
  length <- length(outAve)
  if(length == 100){
  }
  if(length < 100){
    vecNA <- rep(NA, 100-length)
    outAve<- c(vecNA, outAve)
  }
  if(length > 100) {
    outAve <- tail(outAve, 100)
  }
  
  positionAve <- as.data.frame(outAve)
  rm(out, outAve)
  
  newName <- i[1, "X8"]
  colnames(positionAve)<- newName
  
  
  
  if(exists("outFinalLeft")==TRUE){
    outFinalLeft <- cbind(outFinalLeft, positionAve)
  }else{
    outFinalLeft <- positionAve
  }
  
}

#import the right side of the midpoint timing domain region to transition zone

dfRight <- read_tsv("noct15.100kbS.right.bed", col_names = FALSE)

#make a list based on regions
rightSplit <- split(dfRight, dfRight$X8)

#for each region, resize and scale it make it fit a standard sized window.
for(i in rightSplit){ 
  dfw <- i
  lengthRight <- round(length(dfw$X4)/100)
  
  splitVecs<- split(dfw$X4, ceiling(seq_along(dfw$X4)/lengthRight))
  
  for (n in splitVecs){
    meanN<- mean(n)
    if(exists("outAve")==TRUE){
      outAve<- c(outAve, meanN)
    }else{
      outAve <- meanN
    }
  }
  
  
  length <- length(outAve)
  if(length == 100){
  }
  if(length < 100){
    vecNA <- rep(NA, 100-length)
    outAve<- c(outAve, vecNA)
  }
  if(length > 100) {
    outAve <- tail(outAve, 100)
  }
  
  positionAve <- as.data.frame(outAve)
  rm(out, outAve)
  
  newName <- i[1, "X8"]
  colnames(positionAve)<- newName
  
  
  
  if(exists("outFinalRight")==TRUE){
    outFinalRight <- cbind(outFinalRight, positionAve)
  }else{
    outFinalRight <- positionAve
  }
  
}


left <- as.data.frame(t(outFinalLeft))
right <- as.data.frame(t(outFinalRight))

colnames(left)<- c(1:100)
colnames(right)<- c(101:200)


#To prevent outlier positions from skewing to presentation of the data, set the maximul signal to 100. This is based on
#the average signal across all time point equalling ~100. Without setting this feature, repeat regions in the genome continue
# to increase in sequening depth, thus diltuing the signal increase in these regions, as well as consistent increase in the regions
# of interest over the time point prevents appreiating the filling in ans expansion of sigal due to the scale increases as well.

left[left > 100] <- "100"
right[right > 100] <- "100"

left$ID <- row.names(left)
right$ID <- row.names(right)

# Join the left and right halfs of the data, and sort by signal intensity to determine the plto order and put the data into a plotable form.
merged <- full_join(left, right, by = c("ID"="ID"))
row.names(merged)<-merged$ID
mergedMeans <- select(merged, -ID)
mergedMeans <- data.matrix(mergedMeans)

plotOrder <- rowMeans(mergedMeans, na.rm = TRUE)

plotData <- melt(merged, id.vars = "ID")
plotData$variable <- as.integer(plotData$variable)

dfPO<- as.data.frame(plotOrder)
dfPO$ok <- row.names(dfPO)

dfPO <- arrange(dfPO, plotOrder)
dfPO$orderFinal <-row.names(dfPO)


plotData <- full_join(plotData, dfPO, by = c("ID"="ok"))
plotData$orderFinal <- as.integer(plotData$orderFinal)
plotData$variable <- as.integer(plotData$variable)
plotData$value <- as.numeric(plotData$value)
plotData <- arrange(plotData, orderFinal)

#Plot the data to generate the heatmap

p <- ggplot(plotData, aes(y=orderFinal, x=variable)) +
  geom_tile(aes(fill = value), na.rm = FALSE)+
  scale_fill_gradient2(low="white", high="red", limit = c(0, 100), na.value = "#cccccc") +
  theme(axis.text.x = element_text(size = 0),
        axis.text.y = element_text(size = 0, color = "black"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        panel.border = element_rect(color = "black", fill=NA, size=0.5),
        panel.background = element_blank(),
        axis.ticks = element_line(size = 0),
        plot.margin = margin(r=-0.5)) +
  labs(x="ERD or LRD", y="Timing Domain", fill = "")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))



#Import the ERD annotation to mark their location on the final heatmap
erd <- read_tsv("erd.midpoint.bed", col_names = FALSE)
erd$X10 <- 1
erdWorking <- select(erd, X9, X10)
colnames(erdWorking)<- c("X9", "ERD")

plotData2 <- full_join(erdWorking, dfPO, by = c("X9"="ok"))
plotData2$orderFinal <- as.integer(plotData2$orderFinal)
plotDataFinal <- plotData2[complete.cases(plotData2), ]

pERD <- ggplot(plotDataFinal, aes(ERD, orderFinal)) +
  geom_tile(color = "blue") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 10),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(l=-0.5, r=1))+
  scale_x_discrete(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))


library(ggpubr)

#Plot the final figure
finalp<- ggarrange(p,pERD, 
                   ncol = 2, nrow = 1,  align = "v", 
                   widths = c(4, 1), heights = c(1, 1), 
                   legend = "left", 
                   common.legend = TRUE)

