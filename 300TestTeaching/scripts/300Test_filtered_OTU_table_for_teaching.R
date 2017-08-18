# This script does three things:  
# 1) It takes the DAMe filtered OTU table and the OTU representative sequences and merges them into an OTU table that has been filtered by taxonomy (here, Arthropoda prob >=0.80 only).  
# 2) Then it carries out a phyloseq analysis of OTU sizes and helps you choose a minimum OTU size, and filters out the smaller OTUs
# 3) Finally, it filters out the small OTUs from the OTU table (OTU x sample) and also creates an R object for downstream community analysis (sample X OTU). 
# 
# 
library("vegan"); packageVersion("vegan")
library("car")
library("RColorBrewer") # for colours
library("phyloseq"); packageVersion("phyloseq")
library("data.table"); packageVersion("data.table")
library("tidyverse"); packageVersion("tidyverse") # includes all data-formatting packages as one big package (e.g. dplyr, tidyr, ggplot2, readr, readxl, tibble, and others)
# library(dplyr); packageVersion("dplyr")  # for some reason, tidyverse doesn't let me have access to dplyr

# rm(list=ls())

# set your path to the OTU_tables folder that was the product of the DAMe pipeline
path_name <- file.path("/Users/Negorashi2011/Xiaoyangmiseqdata/MiSeq_20170410/300Test/scripts/300TestTeaching/analysis/OTUs_min2PCRs_min4copies_2017-08-17_time-1634/OTU_tables/")
setwd(path_name)
getwd()

taxcolnames <- c("OTU","root","root_tax","root_prob","superkingdom","superkingdom_tax","superkingdom_prob","phylum","phylum_tax","phylum_prob","class","class_tax","class_prob","order","order_tax","order_prob","family","family_tax","family_prob","genus","genus_tax","genus_prob","species","species_tax","species_prob")

# import RDP taxonomies and format. This is the output from RDP Classifier, and it only includes the OTU name (e.g. OTU1) and its Classifier-assigned taxonomies

taxonomies_B_97 <- read.table("table_300test_B_97.RDPmidori_Arthropoda.txt", header=F, sep="\t")
View(taxonomies_B_97)

taxonomies_B_97 <- setNames(taxonomies_B_97, taxcolnames)
View(taxonomies_B_97)


# import DAMe OTU tables. These tables have NOT YET been filtered to include only Arthropoda. 

otutable_B_97 <- read.table("table_300test_B_97.txt", header=T, sep = "\t")

otutable_B_97 <- otutable_B_97 %>% dplyr::mutate(OTUreadtot=Hhmlbody+hlllbody+hlllleg+Hhmlleg+hhhlleg+hhhlbody+mmmmbody+mmmmleg) 
	
otutable_B_97 <- otutable_B_97 %>% dplyr::select(OTU,PC,OTUreadtot,Hhmlbody,hlllbody,hlllleg,Hhmlleg,hhhlleg,hhhlbody,mmmmbody,mmmmleg,Seq) # create OTU table with only these columns
View(otutable_B_97)

# join the two tables using right_join(), column "OTU" is coerced to char. This filters the OTU table to include only the OTUs that are in the files:  table_300test_A97.RDPmidori_Arthropoda.txt. Thus, only Arthropoda-assigned OTUs are kept. 
otutablefull_B_97 <- right_join(otutable_B_97, taxonomies_B_97, by="OTU")
View(otutablefull_B_97)

# remove the orig OTU and taxonomies files
rm(taxonomies_B_97)
rm(otutable_B_97)

# Now you have a taxonomically filtered OTU table.  


#################################
# phyloseq filtering by OTU size. This step removes "small" OTUs, which are probably artefacts of PCR and sequencing error, i.e. echo OTUs, which should have been clustered into the legitimate OTUs. What is "small"?  It is a bit subjective, but the phyloseq method creates a rarefaction curve, showing how many OTUs are at each size class. "Small" OTUs are suggested when there are *very many* OTUs at small sizes. See the graph created below to understand this.
# Read in an OTU table and convert from OTU X Sample to Sample X OTU


communityAll_t <- otutablefull_B_97 %>% dplyr::select(one_of(c("OTU","PC", "Hhmlbody","hlllbody","hlllleg","Hhmlleg","hhhlleg","hhhlbody","mmmmbody","mmmmleg")))  # note that in some OTU tables, the order of variables is OTU, PC, then the samples.

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



#### phyloseq filtering
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
pCumSum + scale_x_continuous(breaks = scales::pretty_breaks(n = 25)) + scale_y_continuous(breaks = scales::pretty_breaks(n = 25)) # the whole curve. not terribly useful.
# pCumSum + scale_x_continuous(breaks = scales::pretty_breaks(n = 50), limits = c(0, 10000)) + scale_y_continuous(breaks = scales::pretty_breaks(n = 25))
# pCumSum + scale_x_continuous(breaks = scales::pretty_breaks(n = 25), limits = c(0, 500)) + scale_y_continuous(breaks = scales::pretty_breaks(n = 25))
pCumSum + scale_x_continuous(breaks = scales::pretty_breaks(n = 25), limits = c(0, 100)) + scale_y_continuous(breaks = scales::pretty_breaks(n = 25))

cat("There are", dim(communityAll)[2], "OTUs before phyloseq filtering. \n")

threshold_otu_size <- 6
# This dataset has basically no artefactual OTUs left.  I will have to show you pictures from other libraries with. 
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

assign(paste0("OTUtablefull_B_97_minOTU_", threshold_otu_size), left_join(communityAll_t, otutablefull_B_97, by="OTU") %>% select(-one_of("V1","V2","V3","V4","V5","V6","V7","V8","PC")))
View(OTUtablefull_B_97_minOTU_6)

names(OTUtablefull_B_97_minOTU_6)

# The filename has the format:  Otutablefull_A_97_minOTU_6. The last number is the minOTU size.

##### Create a list that holds the sample names and the taxonomy and OTU-size-filtered OTU table. 
assign(paste0("Comm_analysis_list_B_97_minOTU_", threshold_otu_size), list(sample_names.df, communityAll))

# The filename has the format:  Comm_analysis_list_A_97_minOTU_12


#### Congratulations, you now have your final datasets for community analysis. 

# otutablefull_B_97:  The original files from DAMe + RDP Classifier, filtered for taxonomy (typically = Arthropoda at prob ≥ 0.80), in OTU X Sample format, plus RDP Classifier taxonomy)

# Otutablefull_B_97_minOTU_6:  The above file now also filtered for min OTU size, and without the positive control (PC) sample, in OTU X Sample format, plus RDP Classifier taxonomy)

# Comm_analysis_list_B_97_minOTU_6:  The community analysis files in a list:  one with the sample names and one with the OTUs, in Sample X OTU format, no taxonomy)  This is the one that you will likely use the most.  The sample names can be joined with sample information into a larger environment table

###########
# write otutablefull_A_97 tables to disk

otutablelist <- ls(pattern = "otutablefull") # gets all filenames with "otutablefull" in the name
otutablelist
for(i in 1:length(otutablelist))
{
  write.csv(x = get(otutablelist[i]), file = paste0(otutablelist[i], "_with_tax.csv"), row.names = FALSE)
}

# write Otutablefull_A_97_minOTU_12 tables to disk

Otutablelist <- ls(pattern = "OTUtablefull") # gets all filenames with "Otutablefull" in the name
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
# Community analysis
####################################################################

# To read in a list for analysis
communitylist <- readRDS("Comm_analysis_list_B_97_minOTU_6.rds")
env <- communitylist[[1]]
community <- communitylist[[2]]


# do NMDS analysis to see basic patterns ####
bodypart <- setNames(as.data.frame(as.vector(c("body", "body", "leg", "leg", "leg", "body", "body", "leg"))), c("bodypart"))
evenness <- setNames(as.data.frame(as.vector(c("Hhml", "hlll", "hlll", "Hhml", "hhhl", "hhhl", "mmmm", "mmmm"))), c("evenness"))
env <- bind_cols(env, bodypart, evenness)

# species richness for all samples
(sprichness <- specnumber(community, groups = env$sample_names, MARGIN = 1))

# NMDS for all samples

community.jmds <- metaMDS(community, distance = "jaccard", trymax = 20, binary=FALSE)
community.jmds <- metaMDS(community, distance = "jaccard", binary = FALSE, previous.best = community.jmds)
stressplot(community.jmds)


# change the order of the levels in env$bodypart, so that the legend ordered up to down like the points
# RUN ONCE, or the bodypart order will change again
levels(env$bodypart)
env$bodypart <- factor(env$bodypart, levels(env$bodypart)[c(2,1)])
levels(env$bodypart)


lib <- community.jmds

	## extract scores
	sites <- scores(lib, display = "sites")
	spps  <- scores(lib, display = "species")
	
	## compute axis ranges
	xlim <- range(sites[,1], spps[,1])
	ylim <- range(sites[,2], spps[,2])
	
	## color palette
	colorvec <- c("#EF8A62", "#67A9CF")  # from:  brewer.pal(3,"RdBu")
	# (colorvec <- c("red2", "mediumblue"))
	## plot
	plot(lib, ylab="", xlab="", xlim=xlim, ylim=ylim, type="n", scaling = 3, main = "community B. point size ~ species richness")
	points(lib, display = "sites", pch=16, cex=sprichness/30, col=colorvec[env$bodypart])
	with(env, legend("topright", legend = levels(bodypart), bty = "n", col=colorvec, pt.cex=2, pch=16))
	cexnum <- 1
	with(env, ordisurf(lib, sprichness, main="", cex=0.5))
	with(env, ordispider(lib, evenness, cex=cexnum, draw="polygon", col=c("black"), alpha=100, kind="se", conf=0.95, label=TRUE, 
											 show.groups=(c("hlll"))))
	with(env, ordispider(lib, evenness, cex=cexnum, draw="polygon", col=c("black"), alpha=100, kind="se", conf=0.95, label=TRUE,
											 show.groups=(c("Hhml"))))
	with(env, ordispider(lib, evenness, cex=cexnum, draw="polygon", col=c("black"), alpha=100, kind="se", conf=0.95, label=TRUE,
											 show.groups=(c("hhhl"))))
	with(env, ordispider(lib, evenness, cex=cexnum, draw="polygon", col=c("black"), alpha=100, kind="se", conf=0.95, label=TRUE,
											 show.groups=(c("mmmm"))))



# color palette webpages
# http://www.fromthebottomoftheheap.net/2012/04/11/customising-vegans-ordination-plots/
# https://github.com/hallamlab/mp_tutorial/wiki/Taxonomic-Analysis




