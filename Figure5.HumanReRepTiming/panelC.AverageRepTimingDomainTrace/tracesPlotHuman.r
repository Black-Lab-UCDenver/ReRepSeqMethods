
#Import libraries
library(readr)
library(dplyr)
library(ggplot2)

#set the working directory
setwd("/path/to/files")

#load input files
df0 <- read_tsv("noct0.100kbS.erd.lrd.averages.traces.txt", col_names = TRUE)
df10 <- read_tsv("noct10.100kbS.erd.lrd.averages.traces.txt", col_names = TRUE)
df15 <- read_tsv("noct15.100kbS.erd.lrd.averages.traces.txt", col_names = TRUE)
df25 <- read_tsv("noct25.100kbS.erd.lrd.averages.traces.txt", col_names = TRUE)
dfh2cc <- read_tsv("H2CC.100kbS.erd.lrd.averages.traces.txt", col_names = TRUE)

#merge the input files into a single file for plotting
plotData <- rbind(df0, df10, df15, df25, dfh2cc)

#split the data in to erds and lrds to make separate plots.
erd <- filter(plotData, domain=="ERD")
lrd <- filter(plotData, domain=="LRD")

#plot the erd traces
pERD<- ggplot(erd, aes(x = order, y = Value, color = sample, group = sample, linetype = sample))+
  geom_smooth(method = "auto", se = FALSE)+
  scale_color_manual(values = c("0hr" = "#00f7ff", "10hr" = "#0400ff", "15hr"= "#8400ff", "25hr"="#ff0026", "H2CC"="black"))+
  scale_linetype_manual(values=c("0hr" = "solid", "10hr" = "solid", "15hr"= "solid", "25hr"="solid", "H2CC"="twodash"))+
  ylim(0, 150)+
  theme(panel.background = element_rect(fill = "white"), 
        legend.background = element_rect(fill = "white"),
        axis.text.x = element_text(size = 0),
        axis.text.y = element_text(size = 12, color = "black"), 
        panel.border = element_rect(color = "black", fill=NA, size=0.5), 
        axis.ticks = element_line(size = 0), 
        legend.key = element_rect(colour = NA, fill = NA))

#plot the lrd traces.
pLRD<- ggplot(lrd, aes(x = order, y = Value, color = sample, group = sample, linetype = sample))+
  geom_smooth(method = "auto", se = FALSE)+
  scale_color_manual(values = c("0hr" = "#00f7ff", "10hr" = "#0400ff", "15hr"= "#8400ff", "25hr"="#ff0026", "H2CC"="black"))+
  scale_linetype_manual(values=c("0hr" = "solid", "10hr" = "solid", "15hr"= "solid", "25hr"="solid", "H2CC"="twodash"))+
  ylim(0, 150)+
  theme(panel.background = element_rect(fill = "white"), 
        legend.background = element_rect(fill = "white"),
        axis.text.x = element_text(size = 0),
        axis.text.y = element_text(size = 12, color = "black"), 
        panel.border = element_rect(color = "black", fill=NA, size=0.5), 
        axis.ticks = element_line(size = 0), 
        legend.key = element_rect(colour = NA, fill = NA))
