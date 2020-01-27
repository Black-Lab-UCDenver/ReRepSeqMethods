
#load necessaryt libraries
library(readr)
library(dplyr)
library(ggplot2)

# set path to working directory
setwd("/path/to/data")


#read input data
df1 <- read_tsv("cdkbypassAllSamples.bed", col_names = FALSE)


#remove begraph coordinates
df2 <- select(df1, -X1, -X2, -X3)
df3 <- as.data.frame(t(df2))

#perform pca
PC<-prcomp(df3)

#ensure the groups are labeled uniformly
Group = c("Gal", "Gal", "Gal", "Raf", "Raf", "Raf")


group <- data.frame(Group)

PCi<-data.frame(PC$x,Group=group$Group)


#plot the data
plot1 = ggplot(PCi,aes(x=PC1,y=PC2,col=Group))+
  geom_point(size=5)+
  scale_color_manual(values = c("Gal" = "red", "Raf"="black"))+ #your colors here
  theme_classic()+ 
  theme(axis.text = element_text(size = 12))+ # changes axis labels
  theme(axis.text=element_text(colour="black"))+
  theme(axis.title = element_text(size = 20)) # change axis titles

p <- p + theme(text = element_text(size = 10)) # this will change all text size 

