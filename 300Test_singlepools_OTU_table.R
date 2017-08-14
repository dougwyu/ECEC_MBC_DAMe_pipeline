# This script is applied to the singlepools
# This script does three things:  
# 1) It takes the DAMe filtered OTU table and the OTU representative sequences and merges them into an OTU table that has been filtered by taxonomy (here, Arthropoda prob >=0.80 only).  
# 2) Then it carries out a phyloseq analysis of OTU sizes and helps you choose a minimum OTU size, and filters out the smaller OTUs
# 3) Finally, it filters out the small OTUs from the OTU table (OTU x sample) and also creates an R object for downstream community analysis (sample X OTU). 
# 
# 
library(vegan); packageVersion("vegan")
library(car)
library(breakaway); packageVersion("breakaway")
library(iNEXT); packageVersion("iNEXT")
library(boral); packageVersion("boral")
library(mvabund); packageVersion("mvabund")
library(PDcalc); packageVersion("PDcalc")
library(tidyverse) # includes all data-formatting packages as one big package (e.g. dplyr, tidyr, ggplot2, readr, readxl, tibble, and others)
library(dplyr); packageVersion("dplyr")  # for some reason, tidyverse doesn't let me have access to dplyr
library("phyloseq"); packageVersion("phyloseq")
library("data.table"); packageVersion("data.table")
library("ggplot2"); packageVersion("ggplot2")

path_name <- file.path("/Users/Negorashi2011/Xiaoyangmiseqdata/MiSeq_20170410/300Test/analysis/singlepools/")
setwd(path_name)
print(path_name)

folder <- c("A", "B", "C", "D", "E", "F")
pools <- c("1", "2", "3")

taxcolnames <- c("OTU","root","root_tax","root_prob","superkingdom","superkingdom_tax","superkingdom_prob","phylum","phylum_tax","phylum_prob","class","class_tax","class_prob","order","order_tax","order_prob","family","family_tax","family_prob","genus","genus_tax","genus_prob","species","species_tax","species_prob")

# import RDP taxonomies and format
# needed to do a lot of fixing of the RDP taxonomies
for(i in folder)
{
	for(j in 1:3)
	{
		for(k in 97:97)
		{
		assign(paste0("taxonomies_", i, j, "_", k), setNames(read.table(paste0("table_300test_", i, j, "_", k, ".RDPmidori_Arthropoda.txt"), header=F, sep = "\t"), taxcolnames))
		# setNames:  sets column names from the taxcolname vector when making (reading) the dataset.
		 assign(paste0("taxonomies_", i, j, "_", k), dplyr::select(get(paste0("taxonomies_", i, j, "_", k)), -contains("_tax")))
		# get:  gets the value of the named object. So if it's a dataframe, get will get the dataframe
		}
	}
}

# OTU tables had TagN-TagN names changed to the sample names in bash before importing here. The changed OTU tables are called:  table_300test_97_A1_samplenames.txt

# import DAMe OTU tables
for(i in folder)
{
	for(j in 1:3)
	{
		for(k in 97:97)
		{
		assign(paste0("otutable_", i, j, "_", k), read.table(paste0("table_300test_", k, "_", i, j, "_samplenames.txt"), header=T, sep = "\t"))
		
		assign(paste0("otutable_", i, j, "_", k), dplyr::mutate(get(paste0("otutable_", i, j, "_", k)), OTUreadtot=Hhmlbody+hlllbody+hlllleg+Hhmlleg+hhhlleg+hhhlbody+mmmmbody+mmmmleg))
		
		assign(paste0("otutable_", i, j, "_", k), dplyr::select(get(paste0("otutable_", i, j, "_", k)), OTU,PC,OTUreadtot,Hhmlbody,hlllbody,hlllleg,Hhmlleg,hhhlleg,hhhlbody,mmmmbody,mmmmleg,xb1,xb2,xb3,Seq))
		}
	}
}


# join the two tables, column "OTU" is coerced to char
for(i in folder)
{
	for(j in 1:3)
	{
		for(k in 97:97)
		{
		assign(paste0("otutablefull_", i, j, "_", k), right_join(get(paste0("otutable_", i, j, "_", k)), get(paste0("taxonomies_", i, j, "_", k)), by="OTU"))
		}
	}
}

# # remove the orig OTU and taxonomies files
rm(list = ls(pattern = "otutable_"))
rm(list = ls(pattern = "taxonomies"))

# Now you have a taxonomically filtered OTU table.  


#################################
# phyloseq filtering by OTU size

# Read in an OTU table and convert from OTU X Sample to Sample X OTU
# eval(parse(text = paste0("otutablefull_", folder, "_", sim))) # to allow dynamic variable names.  It's not terribly useful, however, to automate the following because one needs to look at the intermediate outputs and make decisions.

sim <- 97
folder <- "A"
pool <- 3

communityAll_t <- eval(parse(text = paste0("otutablefull_", folder, pool, "_", sim))) %>% dplyr::select(one_of(c("OTU","PC","OTUreadtot","xb1","xb2","xb3","Hhmlbody","hlllbody","hlllleg","Hhmlleg","hhhlleg","hhhlbody","mmmmbody","mmmmleg")))  # note that in some OTU tables, the order of variables is OTU, PC, then the samples.

# Observe the Positive control OTUs and filter the OTU table.  
communityAll_t <- communityAll_t %>% arrange(desc(PC))
View(communityAll_t)

# In the singlepool samples, large OTUs have some reads in the PC column, and some of the PC OTU reads are in other sample columns.  Cannot keep only the OTUs where PC = 0, because some of the OTUs where PC > 0 are ones with large, real OTUs. 
# Solution:  keep OTUs where PC = 0 OR OTUReadtot > 500.  OTUreadtot does not include PC reads. 500 is a guess, but should include 

communityAll_t <- communityAll_t %>% filter(PC <= 0 | OTUreadtot > 300)
communityAll_t <- communityAll_t %>% select(-c(PC, OTUreadtot, xb1, xb2, xb3)) # remove the PC sample


communityAll <- t(communityAll_t)
rm(communityAll_t)
colvector <- communityAll[1,] # make a vector of the first row, which has the otu names
communityAll <- as.data.frame(communityAll)
colnames(communityAll) <-  colvector # add the otu names to the column names
communityAll <- communityAll[-1,]
sample_names <- rownames(communityAll)
sample_names.df <- data.frame(sample_names)
# convert the columns to numeric from factor
# http://stackoverflow.com/questions/2288485/how-to-convert-a-data-frame-column-to-numeric-type
communityAll <- sapply(communityAll, function(x) as.numeric(as.character(x))) # sapply applies a function to each column, and the function is:  function(x) as.numeric(as.character(x)).  Cannot convert factors to numeric directly. first convert to character, then to numeric
communityAll <- as.data.frame(communityAll) # then convert to df

# Some OTUs might still have 0 total reads. We remove these OTUs
communityAll <- communityAll[ , colSums(communityAll)>0]


# phyloseq code
TotalCounts <- c(colSums(communityAll))

tdt = data.table(OTUs = colnames(communityAll), TotalCounts = colSums(communityAll), OTU = colnames(communityAll))

ggplot(tdt, aes(TotalCounts)) + 
  geom_histogram() + 
  ggtitle("Histogram of Total Counts per OTU")

taxcumsum = tdt[, .N, by = TotalCounts]
setkey(taxcumsum, TotalCounts)
taxcumsum[, CumSum := cumsum(N)]
# Define the plot
pCumSum = ggplot(taxcumsum, aes(TotalCounts, CumSum)) + 
  geom_point() + 
  xlab("Filtering Threshold:  Minimum Read Number per OTU") +
  ylab("OTUs That Would Be Filtered Out") +
  ggtitle("Number of OTUs that would be filtered out at different minimum OTU size thresholds")
pCumSum + scale_x_continuous(breaks = scales::pretty_breaks(n = 25)) + scale_y_continuous(breaks = scales::pretty_breaks(n = 25))
# pCumSum + scale_x_continuous(breaks = scales::pretty_breaks(n = 50), limits = c(0, 10000)) + scale_y_continuous(breaks = scales::pretty_breaks(n = 25))
# pCumSum + scale_x_continuous(breaks = scales::pretty_breaks(n = 25), limits = c(0, 500)) + scale_y_continuous(breaks = scales::pretty_breaks(n = 25))
pCumSum + scale_x_continuous(breaks = scales::pretty_breaks(n = 25), limits = c(0, 100)) + scale_y_continuous(breaks = scales::pretty_breaks(n = 25))

# I look to see if the curve has has an early near-vertical rise, which indicates a large number of small OTUs. In this dataset, there is not much of one, but there is a bit of a near-vertical rise < OTU size threshold = 13, amounting to thousands of OTUs (best seen with x-range limit of 0- 500). This happens to fit with the PC data, but the point is that there isn't a large number of small OTUs, obviously, because we filtered them out via DAMe

# A1: 8, A2: 8, A3: 8
threshold_otu_size <- 8
communityAll <- communityAll[, colSums(communityAll) >= threshold_otu_size]
rowSums(communityAll) # confirm that all samples (rows) still have non-zero OTUs in them.  This isn't a risk with this dataset, but some datasets have samples with very few, small OTUs.  Removing small OTUs will produce samples (rows) that have almost no data.  If so, you'll wan to remove these here

threshold_sample_size <- 2000 # the 2000 is just an example, indicating samples that have at least 2000 reads summed over all OTUs
communityAll <- communityAll[rowSums(communityAll) >= threshold_sample_size, ] 


#### Make a new OTU table, with taxonomic information. Filter out small OTUs from the original otutablefull_A_96 type files. This dataset also no longer has the PC column. 

communityAll_t <- t(communityAll)
communityAll_t <- as.data.frame(communityAll_t)
communityAll_t <- rownames_to_column(communityAll_t)
communityAll_t <- communityAll_t %>% rename(OTU=rowname)

assign(paste0("Otutablefull_", folder, pool, "_", sim, "_minOTU_", threshold_otu_size), left_join(communityAll_t, get(paste0("otutablefull_", folder, pool, "_", sim)), by="OTU") %>% select(-one_of("1","2","3","4","5","6","7","8","PC","OTUreadtot","xb1","xb2","xb3")))

# The filename has the format:  Otutablefull_A1_97_minOTU_8. The last number is the minOTU size.

##### Create a list that holds the sample names and the taxonomy and OTU-size-filtered OTU table. 
assign(paste0("Comm_analysis_list_", folder, pool, "_", sim, "_minOTU_", threshold_otu_size), list(sample_names.df, communityAll))

# The filename has the format:  Comm_analysis_list_A1_96_minOTU_8


###########################################################################
# Congratulations, you now have your final datasets for community analysis.

# otutablefull_A1_97:  The original files from DAMe + RDP Classifier, filtered for taxonomy (typically = Arthropoda at prob ≥ 0.80), in OTU X Sample format, plus RDP Classifier taxonomy)

# Otutablefull_A_97_minOTU_13:  The above file now also filtered for min OTU size, and without the positive control (PC) sample, in OTU X Sample format, plus RDP Classifier taxonomy)

# Comm_analysis_list_A1_96_minOTU_15:  The community analysis files in a list:  one with the sample names and one with the OTUs, in Sample X OTU format, no taxonomy)

###########
# write Otutablefull_A_97_minOTU_13 tables to disk

Otutablelist <- ls(pattern = "Otutablefull")
Otutablelist
for(i in 1:length(Otutablelist))
{
  write.csv(x = get(Otutablelist[i]), file = paste0(Otutablelist[i], "_filtered_with_tax.csv"), row.names = FALSE)
}

# write otutablefull_A_96 tables to disk

otutablelist <- ls(pattern = "otutablefull")
otutablelist
for(i in 1:length(otutablelist))
{
  write.csv(x = get(otutablelist[i]), file = paste0(otutablelist[i], "_with_tax.csv"), row.names = FALSE)
}

# save Comm_analysis_list_A1_97_minOTU_13 lists to disk, as RDS objects. Can be read in again via readRDS.

commlist <- ls(pattern = "Comm_analysis") # gets all filenames with "Comm_analysis" in the name
commlist

for(i in 1:length(commlist))
{
	saveRDS(get(commlist[i]), file = paste0(commlist[i], ".rds"))
}

# To read in the list for analyses
community <- readRDS("Comm_analysis_list_A1_97_minOTU_13")
env <- community[[1]]
otus <- community[[2]]

identical(Comm_analysis_list_A1_97_minOTU_13, Comm_analysis_list_A1_97_minOTU_13_fromdisk) # checks if the RDS object saved to disk and read back in is the same as the one that i saved
