##########################ReadMe############################
# This code will generate the heatmap found in figure 3 panel b. This analysis was performed 
# to look for signal enrichment in erd and lrd timing domains across known ars origins in yeast. 
# These ars sequences are quite small, in the order of 1kb or smaller, which makes the 
# initation of and subsequent spread of dna replication difficult to monitor. To better montitor
# these features of dna replication, we used a 25kb interval from the middle of each ars sequence,
# producing equally sized 50kb genomic coordinates. This larger interval allowed us to visual 
# both the onset of replication (signal in the middle of a given region) as well as the spread of signal.
# 
# We featured the 25min timepoint in this analysis, this data is provided in this directory.
#
#The relevant inputs for this code are
#1) the 25 minute timepoint bed file
#2) the erd key to match erds to confirmed ars sequences.


#Load the necessary libraries
library(readr)
library(dplyr)
library(ggpubr)
library(reshape2)
library(ggplot2)


#set the working directory, where the input files are saved
setwd("/path/to/files")

#Load the 25 min timepoint file
dfIn <- read_tsv("t25.25kars.bg", col_names = FALSE)

#Because this analysis is dependent on resizing and normalizing each ars domain to generate a large matrix, 
# we begin by splitting the input data based on ars.
listSplit <- split(dfIn, dfIn$X8)

# For each ars in the split list, perform the following code. The code that proceedes verifies the regions are the same size
# and either trims the regions or adds the necessary NA values (which show up gray in the figure). The NA values come into play
# when ars sequences are close to the the end sof chromosomes, and thus need "na" pads on either end.
for(i in listSplit){
  
  #make the current ars sequence it's own data frame
  dflist <- as.data.frame(i)
  
  #Get the length of the data frame. Because the bins are all the same size, the length of the list is also equiveleant to the 
  # length of the region in question. In this case, the sample is not bined at all and is at single base pair resolution.
  length <- length(dflist[,1])
 
  #We need to get the starting position of the region. This will come in handy to identify where the region is
  #geographically on the chromosome.
  position <- dflist[1, 6]
  
  #We need the name of the region
  name<- dflist[1, 8]
  
  #If the length of the region is exactly 50kb, which it should be, then select the column with the rerepsignal in it, in preparation
  #to merge this other processed ars regions to make a large matirx
  if(length == 50000){
    
    out<- dflist[,4]
    
  }#If the length of the region is less than 50k, and the region happens to start at the end of the chromosome
  #add NA padding fill in the space on the right side of the region, thus keeping the ars in the middle of the figure.
  if(length < 50000 & position> 1){
      
      vecNA <- rep(NA, 50000-length)
      signal <- dflist[,4]
      out<- c(signal, vecNA)
      
  }# If the length of the region is less than 50k, and the region starts at the begining of the chromosome add
  # NA padding to the left side until the region is 50k in total length.
  if(length < 50000 & position == 1) {

    vecNA <- rep(NA, 50000-length)
    signal <- dflist[,4]
    out<- c(vecNA, signal)
    
  }
  

#Plotting 50k basepairs is not necessary or realistic. To make the data easier to display, 
#this next block of code re-bins the vectors in to 250 bins of equal length and averages the signal over these bins.
#Please note that NA bins are not considered when compouting the average signal for each bin, unless the entire bin is NAs, 
# in which case, the value for that bin is remain NA, and will be displayed in gray in the final figure.
    splitVecs <- split(out, ceiling(seq_along(out)/200))
  
  for (n in splitVecs){
    meanN<- mean(n)
    if(exists("outAve")==TRUE){
      outAve<- c(outAve, meanN)
    }else{
      outAve <- meanN
    }
  }
    
# Now that the ars region has been adjusted to a length of 50k and rebinned to 250, its time to make the final data matrix
# by iteratively adding each region as it's own column in the final matrix as it finishes processing.
  if(exists("outFinal")==TRUE){
    outdf<- as.data.frame(outAve)
    colnames(outdf)<- name
    outFinal<- cbind(outFinal, outdf)
  }else{
    outFinal <- as.data.frame(outAve)
    colnames(outFinal)<- name
  }

  #remove these varibales after each ars region is processed to ensure they are generated based on the following region to be processed.  
  rm(meanN)  
  rm(outAve)
  
}


#The final matrix has been assembled, its time to order the martrix by signal intensity. The following lines of code accomlish this.

# We need to make the ars names their own id column.
outFinal$position <- row.names(outFinal)

#Make all the data by poisition it's own integer. This ensures the data sorts nnumerically instead of alphabetically. 
outFinal$position <- as.integer(outFinal$position)

# To plot the data, it needs to be in "long" data formate.
plotData <- melt(outFinal, id.vars = "position")

# Next we want to sort this data by position. In the step after this, we use the position as the key to add the signal-sorted plot order.
plotData<- arrange(plotData, position)




#Now that we have the data to be plotted in the correct format and we have a key to which we can add additonal variales, its time to
#compute the average signal across the bins to determime the signal-ranked plotting order.

#We don't need the position variable that we created a moment ago, because this is now an interger it will contribute to the calculated average.
outFinal <- select(outFinal, -position)

#Generate the mean values for each region, and don't consider NA values when performing this calculation.
#Remember the NA values occure when the ends of the chromosomes, we dont want to count empty space.
orderDf <- colMeans(outFinal, na.rm = TRUE)
orderDF <- as.data.frame(orderDf)
colnames(orderDF)<-c("arsVal")
orderDF$ARS <- row.names(orderDF)

#Order the data based on it's signal. 
orderDF2 <- arrange(orderDF, arsVal)
orderDF2$plotOrder<- row.names(orderDF2)

#Merge the correclty formated data for plotting with the plot order we just defined.
plotData2 <- full_join(orderDF2, plotData, by = c("ARS"="variable"))
plotData2$plotOrder <- as.integer(plotData2$plotOrder)

#Based on this experient and the timepoint we are using, we expect the early regions to increase in signal first.
#Import the ERD regon file so we can annotate the final heatmap with the ERD position within the ranked list.
erd <- read_tsv("arsConfirmedERDOverlap.bed", col_names = FALSE)
erdWorking <- select(erd, X4, X8)

#Merge ERD list with the data.
erdOrder <- full_join(orderDF2, erdWorking, by = c("ARS" = "X4"))

#Indiacte the ERDs with a "1", so we can set the color for the ERD anootation.
erdOrder$ERD <- ifelse(erdOrder$X8>0, "1", erdOrder$X8)
erdOrder$plotOrder <- as.integer(erdOrder$plotOrder)
erdOrder2 <- filter(erdOrder, ERD>0)


#Make the erd annotation plot
pERD <- ggplot(erdOrder2, aes(ERD, plotOrder)) +
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


#Make the data plot
p <- ggplot(plotData2, aes(position, plotOrder)) +
  geom_tile(aes(fill = value))+
  scale_fill_gradient2(low="white", high="red", limit = c(0, 15), na.value = "#cccccc") +
  theme(axis.text.x = element_text(size = 0),
        axis.text.y = element_text(size = 0, color = "black"), 
        axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10), 
        panel.border = element_rect(color = "black", fill=NA, size=0.5), 
        panel.background = element_blank(), 
        axis.ticks = element_line(size = 0), 
        plot.margin = margin(r=-0.5)) +
  labs(x="ARS +/- 25kb", y="ARS", fill = "")+
  scale_x_discrete(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0)) 




#Position the final data plot and erd plot into the same image, and perfectly align their axises. 
finalp<- ggarrange(p,pERD, 
          ncol = 2, nrow = 1,  align = "v", 
          widths = c(4, 1), heights = c(1, 1), 
          legend = "left", 
          common.legend = TRUE)



  