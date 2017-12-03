## ------------------------------------------------------------------------
library(tidyverse) # includes all data-formatting packages as one big package (e.g. dplyr, tidyr, ggplot2, readr, readxl, tibble, and others)
# library("devtools") 
# install_github("tobiasgf/lulu", force = TRUE)
library(lulu)
# sessionInfo()

## ------------------------------------------------------------------------
# rm(list=ls())
# path_name <- file.path("/Users/Negorashi2011/Xiaoyangmiseqdata/MiSeq_20170410/300Test/analysis/OTUs_min2PCRs_min4copies_2017-11-26_time-1727/OTU_tables/")
# path_name <- file.path("/Users/Negorashi2011/Xiaoyangmiseqdata/MiSeq_20170410/300Test/analysis/OTUs_min1PCRs_min1copies_2017-11-30_time-2148/OTU_tables/")


# set some general options
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE) # if you want to pass arguments to R from the shell (bash) command line
print(args)

otutablefolder <- args[1]

path_name <- file.path(paste0("/Users/Negorashi2011/Xiaoyangmiseqdata/MiSeq_20170410/300Test/analysis/", otutablefolder, "/OTU_tables/"))

setwd(path_name)
# print(path_name)

## ----lulu----------------------------------------------------------------

experiment <- c("A", "B", "C", "D", "E", "F")

for(i in experiment) # input the matchlist files
{
	for(j in 97:97)
	{
		assign(paste0("matchlist_", i, "_", j), read.table(paste0("match_list_", i, j, ".txt"), header = FALSE, as.is = TRUE, stringsAsFactors = FALSE))
	}
}

for(i in experiment) # input the sumaclust otu tables
{
	for(j in 97:97)
	{
		assign(paste0("table_300test_", i, "_", j), read.table(paste0("table_300test_", i, "_", j,".txt"), header = TRUE, as.is = TRUE, stringsAsFactors = FALSE))
	}
}

# names(table_300test_A_97)

for(i in experiment) # select only the columns with OTU read numbers
{
	for(j in 97:97)
	{
	    assign(paste0("table_300test_", i, "_", j, "_lulu"), dplyr::select(get(paste0("table_300test_", i, "_", j)), -Seq)) # create OTU table without the Seq
	}
}

# move OTU names to rownames
for(i in experiment) # select only the columns with OTU read numbers
{
	for(j in 97:97)
	{
	    assign(paste0("table_300test_", i, "_", j, "_lulu"), get(paste0("table_300test_", i, "_", j, "_lulu")) %>% column_to_rownames(var = "OTU")) 
	}
}

## ------------------------------------------------------------------------
# experiment <- c("C") # if i want to test individual datasets
# experiment <- c("A", "B", "C", "D", "E", "F") # to restore orig value of experiment

for(i in experiment) # select only the columns with OTU read numbers
{
	for(j in 97:97)
	{
	    assign("otutable_lulu", get(paste0("table_300test_", i, "_", j, "_lulu")))
	    assign("matchlist_lulu", get(paste0("matchlist_", i, "_", j)))
	    
	    curated_result <- lulu(otutable_lulu, matchlist_lulu)
	    curated_table <- rownames_to_column(curated_result$curated_table, var = "OTU")

	    write.table(curated_table, file = paste0("table_300test_", i, "_", j, "_lulu.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
	    
	    # write.table(curated_result$curated_otus, file = paste0("otus_300test_", i, "_", j, "_lulu.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
	    # write.table(curated_result$original_table, file = paste0("table_300test_", i, "_", j, ".txt"), sep = "\t")
	    rm(otutable_lulu)
	    rm(matchlist_lulu)
	}
}


# default "curated_result <- lulu(otutab, matchlist, minimum_ratio_type = "min", minimum_ratio = 1, minimum_match = 84, minimum_relative_cooccurence = 0.95)"

# curated_result$curated_table
# curated_result$original_table
# curated_result$curated_count
# rm(otutable_lulu)
# rm(matchlist_lulu)

## ----eval=FALSE, include=FALSE-------------------------------------------
## library(knitr)
## setwd("~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/scripts")
## purl("LULU.Rmd")

