
################Read Me#############

####################################
# This is the normalization sciprt used to normalize rerepseq samples by rpm and mitochondria DNA.
# Three inputs are necessary to run this script:
# 1) Bedgraph file, generated by bedtools.
# 2) Output from samtools idxstats for the sample to be normalized
# 3) Output from samtools idxstats for the blacklisted regions for the sample to be normalized.
# 
# This algorithm calculates the RPM based on reads that are not aligned to mitochondrial DNA or black listed regions, 
# and then scales the resulting value based on mitochondria read count.
####################################


# load the necessary libraries
library(readr)
library(dplyr)

#set the working directory where the input files are located
setwd("/path/to/files")

#load the bedgraph an alignment files. dfIn is the bedgraph, 
#dfCounts is the idxstats samtools output for the sample to be normalized, 
#and dfBL is the idxstats samtools output for the blacklisted regions for the sample to be normalized. 
dfIn <- read_tsv("sample.bg", col_names = FALSE) 
dfCounts <- read_tsv("sampleStats.txt", col_names = FALSE)
dfBL <- read_tsv("BlStats.txt", col_names = FALSE)

# We are only working with the main chromosomes, so we want to remove the stats from alternative chromosomes before we calculate our normalization factors.
# If you are apply this to yeast, make sure to change the chromosomes below to match yeast notation.
# Note that the filter function used is selecting for the chromosomes we want to keep, which includes chrM and chrX
dfCountsMain <- filter(dfCounts, X1=="chrM" | X1 == "chr1" | X1 == "chr2" | X1 == "chr3" | X1 == "chr4" | X1 == "chr5" | X1 == "chr6" | X1 == "chr7" | X1 == "chr8" | X1 == "chr9" | X1 == "chr10" | X1 == "chr11" | X1 == "chr12" | X1 == "chr13" | X1 == "chr14" | X1 == "chr15" | X1 == "chr16" | X1 == "chr17" | X1 == "chr18" | X1 == "chr19" | X1 == "chr20" | X1 == "chr21" | X1 == "chr22" | X1=="chrX")

# Once the main chromosomes have been selected, split this file such that chrM is in it's own dataframe.
dfCountF <- filter(dfCountsMain, X1 !="chrM")  
chrM <- filter(dfCounts, X1 =="chrM")

# Sum the total number of reads aligned to the main chromosomes.
totalCount <- sum(dfCountF$X3)

# Sum the total number of reads aligned to chrM.
totalM <- sum(chrM$X3)

# Sum the total number of reads aligned to blacklisted regions.
totalBL <- sum(dfBL$X3)


# Perform the RPM calculation. Note we are not including chrM in this calculation, and the black listed reads are also removed.
dfIn$X5 <- dfIn$X4/((totalCount-totalBL)/1000000)

# Perform the mitochondrial DNA scaling function. 
dfIn$X6 <- dfIn$X5/(totalM/(totalCount-totalBL))

# Filter the bedgraph file to select for the main chromosomes. 
dfInMain <- filter(dfIn, X1=="chrM" | X1 == "chr1" | X1 == "chr2" | X1 == "chr3" | X1 == "chr4" | X1 == "chr5" | X1 == "chr6" | X1 == "chr7" | X1 == "chr8" | X1 == "chr9" | X1 == "chr10" | X1 == "chr11" | X1 == "chr12" | X1 == "chr13" | X1 == "chr14" | X1 == "chr15" | X1 == "chr16" | X1 == "chr17" | X1 == "chr18" | X1 == "chr19" | X1 == "chr20" | X1 == "chr21" | X1 == "chr22" | X1=="chrX")

# Now, select the columns containig the bedgraph coordinates and the final normalized data.
dfOut <- select(dfInMain, X1, X2, X3, X6)

# Regions with 0 read coverage will often be listed as a -infinity, instead of 0, when applying the above the calculations. Make sure these
# are correctly reported as 0 using the following ifelse statement.
dfOut$X7 <- ifelse(dfOut$X6<0, "0", dfOut$X6)

# The outout of the ifelse statement above turns the numbers into characters, use the following script to turn these back into numbers.
dfOutFinal <- select(dfOut, X1, X2, X3, X7)
dfOutFinal$X7 <- as.numeric(as.character(dfOutFinal$X7))

# Save the file.
write_tsv(dfOutFinal, "Sample.rpmMito.bedgraph", col_names = FALSE)