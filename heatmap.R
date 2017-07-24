# This loop runs through all the folders and pools within folders to calculate a heatmap.  


library(dplyr)
library(tidyr)
folder <- c("A", "B", "C", "D", "E", "F")

for(i in folder)
  for(j in 1:3)
{
  folder <- i
  poolnumber <- j
  path_name <- file.path("/Users/Negorashi2011/Xiaoyangmiseqdata/MiSeq_20170410/300Test/data/seqs", paste("folder", folder, sep=""), paste("pool", poolnumber, sep=""))
  
  setwd(path_name)
  print(path_name)
  
  pool <- read.table("SummaryCounts.txt", header=F, sep = "\t")
  pool <- pool %>% select(-V3) # delete the unique sequences column
  poolspread <- pool %>% spread(V2, V4) # tidyr::spread
  
  row.names(poolspread) <- poolspread$V1
  poolspread <- poolspread %>% select(-V1)
  poolspread <- data.matrix(poolspread)
  
  pdf(paste("heatmap_", "Folder", folder, "_Pool", poolnumber, ".pdf", sep=""))
  poolspread_heatmap <- heatmap(poolspread, Rowv=NA, Colv=NA, revC = TRUE, col = cm.colors(256), scale="column", margins=c(5,10))
  dev.off()
}



# https://flowingdata.com/2010/01/21/how-to-make-a-heatmap-a-quick-and-easy-solution/



