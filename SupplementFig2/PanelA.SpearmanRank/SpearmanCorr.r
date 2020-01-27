library(readr)
library(dplyr)
library(ggplot2)

setwd("/path/to/data")


df1 <- read_tsv("allSamplesUnion.bedgraph", col_names = FALSE)


df2 <- select(df1, -X1, -X2, -X3)
df3 <- as.data.frame(t(df2))


colnames(df2) = c("s1.0hr", "s2.0hr", "s3.0hr", "s1.10hr", "s2.10hr", "s3.10hr", "s1.15hr", "s2.15hr", "s3.15hr", "s1.25hr", "s2.25hr", "s3.25hr")
subSet <- select(df2, s1.25hr, s2.25hr, s3.25hr)
m<- cor(subSet, method="spearman")

m<- as.data.frame(m)

write_tsv(m, "humanNocSpearmanRankReplicates.txt", col_names = TRUE)

