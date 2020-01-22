


##########Read Me####################
#####################################
# This script is used to smooth rerepseq data, and requires the RPM normalized and mitochondria DNA scalled, re-binned begraph.
# The smoothing window function relies on the zoo package, and this will smooth the data based on the size of the window. This 
# means if the bedgraph input files was binned to 1000bp for example, and a 100 row smoothing window is used, this becomes a 100kb
# smoothing window. Please double check the width of your bins is correct, and that the width of the smoothing window below is also 
# the correct width.



# Load the necessary libraries
library(readr)
library(dplyr)
library(zoo)


# Set the working directory to the location of your files.
setwd("/path/to/files")

# Read in the standardized binned bedgraph
dfIn <- read_tsv("inSample.bedgraph", col_names=FALSE)

# Each chromosome needs to be smoothed indivdually, this function splits the sample into a list of files, 
# with each file being an individual chromosome
dfsplit <- split(dfIn, dfIn$X1)

# This function applies the smoothing window to each chromosome in the list.
for(i in dfsplit){
  
  # Ensure the data is in a dataframe format
  inData <- as.data.frame(i)
  # Ensure the bins are in order
  inData<- arrange(inData, X2)
  # Generate a new column containing the rolling mean (aka smoothed data).
  # Please make sure the window size is correct. By aligning the results to center, the ends of the choromoses are symetrically padded. 
  inData$X5 <- rollmean(inData$X4, 100, align = "center", na.pad = TRUE)
  
  #This line of code extracts the chromosome name, renames the dataframe just processed to the name of the chromosome, 
  #and saves this in the environment.
  nextData <- as.data.frame(inData)
  chrName <- nextData[1,1]
  assign(paste(chrName), nextData)
  
  
}

 # Now that all the chromosomes are processed and smoothed, reassembled the full bedgraph.
 # Because the chromosomes are named, you can simply use the rbind function to compile them in any order desired.
 # Note that the chromosome names will beed to be changed to reflect the speices being processed.
out <- rbind(chrI, chrII, chrIII, chrIV, chrV, chrVI, chrVII, chrVIII, chrIX, 
             chrX, chrXI, chrXII, chrXIII, chrXIV, chrXV, chrXVI)

# Select for complete cases only.
outFinal <- out[complete.cases(out), ]

# Select the columns that correspond to the bedgraph coordinates and the smootheed data.
outFinal <- select(outFinal, X1, X2, X3, X5)

# Save the final file.
write_tsv(outFinal, "SampleOut.bg", col_names = FALSE) 
