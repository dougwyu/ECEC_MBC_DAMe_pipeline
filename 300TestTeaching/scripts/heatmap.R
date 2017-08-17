# Generates a heatmap from a SummaryCounts.txt file 

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)


library(dplyr)
library(tidyr)

  pool <- read.table("SummaryCounts.txt", header=F, sep = "\t")
  pool <- pool %>% select(-V3) # delete the unique sequences column
  poolspread <- pool %>% tidyr::spread(V2, V4) 
  
  row.names(poolspread) <- poolspread$V1
  poolspread <- poolspread %>% select(-V1)
  poolspread <- data.matrix(poolspread)
  
  pdf("heatmap.pdf")
  poolspread_heatmap <- heatmap(poolspread, Rowv=NA, Colv=NA, revC = TRUE, col = cm.colors(256), scale="column", margins=c(5,10))
  dev.off()


# https://flowingdata.com/2010/01/21/how-to-make-a-heatmap-a-quick-and-easy-solution/



