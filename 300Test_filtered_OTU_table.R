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
library("RColorBrewer") # for colours
# library("devtools") 
# devtools::install_github("grunwaldlab/metacoder")
library(metacoder)
sessionInfo()

# rm(list=ls())
path_name <- file.path("/Users/Negorashi2011/Xiaoyangmiseqdata/MiSeq_20170410/300Test/analysis/OTUs_min2PCRs_min4copies_2017-08-09_time-1249/OTU_tables/")
setwd(path_name)
print(path_name)

folder <- c("A", "B", "C", "D", "E", "F")

taxcolnames <- c("OTU","root","root_tax","root_prob","superkingdom","superkingdom_tax","superkingdom_prob","phylum","phylum_tax","phylum_prob","class","class_tax","class_prob","order","order_tax","order_prob","family","family_tax","family_prob","genus","genus_tax","genus_prob","species","species_tax","species_prob")

# import RDP taxonomies and format. This is the output from RDP Classifier, and it only includes the OTU name (e.g. OTU1) and its Classifier-assigned taxonomies. 
for(i in folder)
{
	for(j in 96:97)
	{
		assign(paste0("taxonomies_", i, "_", j), setNames(read.table(paste0("table_300test_", i, "_", j, ".RDPmidori_Arthropoda.txt"), header=F, sep = "\t"), taxcolnames))
		# setNames:  sets column names from the taxcolname vector when making (reading) the dataset.
		 assign(paste0("taxonomies_", i, "_", j), dplyr::select(get(paste0("taxonomies_", i, "_", j)), -contains("_tax")))
		# get:  gets the value of the named object. So if it's a dataframe, get will get the dataframe
	}
}

# import DAMe OTU tables. These tables have NOT YET been filtered to include only Arthropoda. 
for(i in folder)
{
	for(j in 96:97)
	{
		assign(paste0("otutable_", i, "_", j), read.table(paste0("table_300test_", i, "_", j, ".txt"), header=T, sep = "\t"))
		
		assign(paste0("otutable_", i, "_", j), dplyr::mutate(get(paste0("otutable_", i, "_", j)), OTUreadtot=Hhmlbody+hlllbody+hlllleg+Hhmlleg+hhhlleg+hhhlbody+mmmmbody+mmmmleg)) # create new column to sum read numbers over all soups.
		
		assign(paste0("otutable_", i, "_", j), dplyr::select(get(paste0("otutable_", i, "_", j)), OTU,Hhmlbody,hlllbody,hlllleg,Hhmlleg,hhhlleg,hhhlbody,mmmmbody,mmmmleg,PC,OTUreadtot,Seq)) # create OTU table with only these columns
	}
}


# join the two tables using right_join(), column "OTU" is coerced to char. This filters the OTU table to include only the OTUs that are in the files:  table_300test_A97.RDPmidori_Arthropoda.txt. Thus, only Arthropoda-assigned OTUs are kept. 
for(i in folder)
{
	for(j in 96:97)
	{
		assign(paste0("otutablefull_", i, "_", j), right_join(get(paste0("otutable_", i, "_", j)), get(paste0("taxonomies_", i, "_", j)), by="OTU"))
	}
}

# # remove the orig OTU and taxonomies files
rm(list = ls(pattern = "otutable_"))
rm(list = ls(pattern = "taxonomies"))

# Now you have taxonomically filtered OTU tables.  


#################################
# phyloseq filtering by OTU size. This step removes "small" OTUs, which are probably artefacts of PCR and sequencing error, i.e. echo OTUs, which should have been clustered into the legitimate OTUs. What is "small"?  It is a bit subjective, but the phyloseq method creates a rarefaction curve, showing how many OTUs are at each size class. "Small" OTUs are suggested when there are *very many* OTUs at small sizes. See the graphs created below to understand this.

# Read in an OTU table and convert from OTU X Sample to Sample X OTU
# eval(parse(text = paste0("otutablefull_", folder, "_", sim))) # to allow dynamic variable names.  It's not useful to make a loop of the following code because one needs to look at the intermediate outputs and make decisions. Luckily, there are only 6 (A-F) to run, per sumaclust similarity percentage (e.g. 97%)

folder <- "F"
sim <- 96

cat("Processing sample ", folder,"_",sim, "\n", sep="")

communityAll_t <- get(paste0("otutablefull_", folder, "_", sim)) %>% dplyr::select(one_of(c("OTU","PC", "Hhmlbody","hlllbody","hlllleg","Hhmlleg","hhhlleg","hhhlbody","mmmmbody","mmmmleg")))  # note that in some OTU tables, the order of variables is OTU, PC, then the samples.

# Observe the Positive control OTUs and filter the OTU table.  
communityAll_t <- communityAll_t %>% arrange(desc(PC))
View(communityAll_t)

# If, for instance, the PC OTUs are only present in the PC sample, then filter out those OTUs. Also, the number of reads per true PC OTU (there are 4 in 300Test) versus the false, smaller, "echo OTUs" gives us guidance on what an artefactual OTU size is (although we generally expect "echo OTUs" in the PC sample to be bigger than echo OTUs in real samples, because there are fewer species in PC samples). 
communityAll_t <- communityAll_t %>% filter(PC <= 0)
communityAll_t <- communityAll_t %>% select(-PC) # remove the PC sample
View(communityAll_t)

#### Now transpose to make canonical OTU tables (sample X OTU format) for community analysis
communityAll <- t(communityAll_t)
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
# Some OTUs might have 0 total reads. We remove these OTUs if they exist.
communityAll <- communityAll[ , colSums(communityAll)>0]
View(communityAll)
rm(communityAll_t)

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
  ylab("Number of OTUs That Would Be Filtered Out") +
  ggtitle("Number of OTUs that would be filtered out at different minimum OTU read-numbers")
# pCumSum + scale_x_continuous(breaks = scales::pretty_breaks(n = 25)) + scale_y_continuous(breaks = scales::pretty_breaks(n = 25)) # the whole curve. not terribly useful.
# pCumSum + scale_x_continuous(breaks = scales::pretty_breaks(n = 50), limits = c(0, 10000)) + scale_y_continuous(breaks = scales::pretty_breaks(n = 25))
# pCumSum + scale_x_continuous(breaks = scales::pretty_breaks(n = 25), limits = c(0, 500)) + scale_y_continuous(breaks = scales::pretty_breaks(n = 25))
pCumSum + scale_x_continuous(breaks = scales::pretty_breaks(n = 25), limits = c(0, 100)) + scale_y_continuous(breaks = scales::pretty_breaks(n = 25))

cat("There are", dim(communityAll)[2], "OTUs before phyloseq filtering. \n")

threshold_otu_size <- 6
# I look to see if the curve has has an early near-vertical rise, which indicates a large number of small OTUs. For example, in the A_97 dataset, there is not much of one, but there is a bit of a near-vertical rise < OTU size threshold = 12, amounting to ~10 OTUs that would be filtered out (best seen with x-range limit of 0-100). This happens to fit with the PC data, but the point is that there isn't a large number of small OTUs, obviously, **because we already filtered them out via DAMe**. This exercise is subjective, but as a general approach, i choose a threshold number equal to the tangent's x-intercept. This gets rid of the OTUs below the steepest climb in OTU numbers. 
communityAll <- communityAll[, colSums(communityAll) >= threshold_otu_size]
rowSums(communityAll) # confirm that all samples (rows) still have non-zero OTUs in them.  This isn't a risk with this dataset, but some datasets have samples with very few, small OTUs.  Removing small OTUs will produce samples (rows) that have almost no data.  If so, you'll wan to remove these here


#### Create a new full OTU table (OTU X sample), including the taxonomic information and filtering out the phyloseq-determined small OTUs from the original otutablefull_A1_97 type files (using left_join()). This dataset also no longer has the PC column. 
communityAll_t <- t(communityAll)
communityAll_t <- as.data.frame(communityAll_t)
communityAll_t <- rownames_to_column(communityAll_t)
communityAll_t <- communityAll_t %>% rename(OTU=rowname)
View(communityAll_t)
names(communityAll_t)

assign(paste0("Otutablefull_", folder, "_", sim, "_minOTU_", threshold_otu_size), left_join(communityAll_t, get(paste0("otutablefull_", folder, "_", sim)), by="OTU") %>% select(-one_of("V1","V2","V3","V4","V5","V6","V7","V8","PC")))

names(get(paste0("Otutablefull_", folder, "_", sim, "_minOTU_", threshold_otu_size)))
# The filename has the format:  Otutablefull_A_97_minOTU_12. The last number is the minOTU size.

##### Create a list that holds the sample names and the taxonomy and OTU-size-filtered OTU table. 
assign(paste0("Comm_analysis_list_", folder, "_", sim, "_minOTU_", threshold_otu_size), list(sample_names.df, communityAll))

# The filename has the format:  Comm_analysis_list_A_97_minOTU_12

########################################################################
###### Now go back up and redo from folder <- "A" to do all six folders
########################################################################



# Congratulations, you now have your final datasets for community analysis. 

# otutablefull_A_96:  The original files from DAMe + RDP Classifier, filtered for taxonomy (typically = Arthropoda at prob ≥ 0.80), in OTU X Sample format, plus RDP Classifier taxonomy)

# Otutablefull_A_96_minOTU_15:  The above file now also filtered for min OTU size, and without the positive control (PC) sample, in OTU X Sample format, plus RDP Classifier taxonomy)

# Comm_analysis_list_A_96_minOTU_15:  The community analysis files in a list:  one with the sample names and one with the OTUs, in Sample X OTU format, no taxonomy)  This is the one that you will likely use the most.  The sample names can be joined with sample information into a larger environment table

###########
# write otutablefull_A_97 tables to disk

otutablelist <- ls(pattern = "otutablefull") # gets all filenames with "otutablefull" in the name
otutablelist
for(i in 1:length(otutablelist))
{
  write.csv(x = get(otutablelist[i]), file = paste0(otutablelist[i], "_with_tax.csv"), row.names = FALSE)
}

# write Otutablefull_A_97_minOTU_12 tables to disk

Otutablelist <- ls(pattern = "Otutablefull") # gets all filenames with "Otutablefull" in the name
Otutablelist
for(i in 1:length(Otutablelist))
{
  write.csv(x = get(Otutablelist[i]), file = paste0(Otutablelist[i], "_filtered_with_tax.csv"), row.names = FALSE)
}


# save Comm_analysis_list_A_97_minOTU_12 lists to disk, as RDS objects. Can be read in again via readRDS.

commlist <- ls(pattern = "Comm_analysis") # gets all filenames with "Comm_analysis" in the name
commlist
for(i in 1:length(commlist))
{
	saveRDS(get(commlist[i]), file = paste0(commlist[i], ".rds"))
}



####################################################################

# To read in the list for analyses
community <- readRDS("Comm_analysis_list_A_96_minOTU_15.rds")
env <- community[[1]]
otus <- community[[2]]

# identical(Comm_analysis_list_A_97_minOTU_12, Comm_analysis_list_A_97_minOTU_12_restored) # identical() checks if the RDS object saved to disk and read back in is the same as the one that was saved. 
