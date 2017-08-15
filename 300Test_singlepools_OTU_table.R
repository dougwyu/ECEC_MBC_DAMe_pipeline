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

# Now you have taxonomically filtered OTU tables.  


#################################
# phyloseq filtering by OTU size

# Read in an OTU table and convert from OTU X Sample to Sample X OTU
# eval(parse(text = paste0("otutablefull_", folder, "_", sim))) # to allow dynamic variable names.  It's not terribly useful, however, to automate the following because one needs to look at the intermediate outputs and make decisions.

sim <- 97
folder <- "A"
pool <- 2

communityAll_t <- get(paste0("otutablefull_", folder, pool, "_", sim)) %>% dplyr::select(one_of(c("OTU","PC","xb1","xb2","xb3","OTUreadtot","Hhmlbody","hlllbody","hlllleg","Hhmlleg","hhhlleg","hhhlbody","mmmmbody","mmmmleg")))  # note that in some OTU tables, the order of variables is OTU, PC, then the samples.

# this used the eval(parse(text))) syntax, which is not as good as get()
# communityAll_t <- eval(parse(text = paste0("otutablefull_", folder, pool, "_", sim))) %>% dplyr::select(one_of(c("OTU","PC","xb1","xb2","xb3","OTUreadtot","Hhmlbody","hlllbody","hlllleg","Hhmlleg","hhhlleg","hhhlbody","mmmmbody","mmmmleg")))  # note that in some OTU tables, the order of variables is OTU, PC, then the samples.

# Observe the Positive control OTUs and filter the OTU table.  
communityAll_t <- communityAll_t %>% arrange(desc(PC))
View(communityAll_t)

# In B1_97, OTUs 3, 17, 20, and 59 have heavy contamination in xb1, so i remove these too
# communityAll_t <- communityAll_t %>% filter(OTU != "OTU3" & OTU != "OTU17" & OTU != "OTU20" & OTU != "OTU59")
# In D1_97, OTUs 73, 72, and 117 are almost certainly PC OTUs that also show up in some of the samples, so need to remove specifically
# communityAll_t <- communityAll_t %>% filter(OTU != "OTU73" & OTU != "OTU72" & OTU != "OTU117")
# In E1_97, sample mmmmleg had very low numbers of reads compared to the other samples. 
# View(communityAll_t)


# In the singlepool samples, large OTUs have some reads in the PC column, and some of the PC OTU reads are in other sample columns.  Cannot keep only the OTUs where PC = 0, because some of the OTUs where PC > 0 are ones with large, real OTUs. 
# Solution:  keep OTUs where PC = 0 OR OTUReadtot > 50.  OTUreadtot does not include PC reads. 50 is a guess, but is set to include OTUs that show up with at least a few dozen reads

communityAll_t <- communityAll_t %>% filter(PC <= 0 | OTUreadtot > 50)
View(communityAll_t)

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
# pCumSum + scale_x_continuous(breaks = scales::pretty_breaks(n = 25)) + scale_y_continuous(breaks = scales::pretty_breaks(n = 25))
# pCumSum + scale_x_continuous(breaks = scales::pretty_breaks(n = 50), limits = c(0, 10000)) + scale_y_continuous(breaks = scales::pretty_breaks(n = 25))
# pCumSum + scale_x_continuous(breaks = scales::pretty_breaks(n = 25), limits = c(0, 500)) + scale_y_continuous(breaks = scales::pretty_breaks(n = 25))
pCumSum + scale_x_continuous(breaks = scales::pretty_breaks(n = 25), limits = c(0, 100)) + scale_y_continuous(breaks = scales::pretty_breaks(n = 25))

# I look to see if the curve has has an early near-vertical rise, which indicates a large number of very small OTUs. I set the min OTU size threshold to be the x-intercept of the tangent to the curve in the limits=c(0,100) plot. In the singlepool datasets, there is a large number of small OTUs, and my choice of threshold is 8 for most sample pools.

threshold_otu_size <- 8
communityAll <- communityAll[, colSums(communityAll) >= threshold_otu_size]

rowSums(communityAll) # Check that all samples (rows) still have sufficient numbers of reads in them.  This isn't a risk with this dataset, but some datasets have samples with very few, small OTUs.  Removing small OTUs will sometimes produce samples (rows) that have almost no data, because that sample failed during PCR or DNA extraction.  

# In sample A2, we find that sample 7 (mmmmbody) has very few total reads (77 total). This sample was mistakenly not included during PCR, so even these 77 reads represent Illumina tag jump or cross-talk contamination.

#### Make a new full OTU table, including the taxonomic information and filtering out the small OTUs from the original otutablefull_A1_97 type files. This dataset also no longer has the PC column. 
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

# Otutablefull_A1_97_minOTU_13:  The above file but now ALSO filtered for min OTU size using the phyloseq method, and without the positive control (PC) sample, in OTU X Sample format, plus RDP Classifier taxonomy)

# Comm_analysis_list_A1_97_minOTU_15:  The community analysis files in a list:  one with the sample names and one with the OTUs, in Sample X OTU format, no taxonomy)

###########
# write otutablefull_A_96 tables to disk

otutablelist <- ls(pattern = "otutablefull")
otutablelist
for(i in 1:length(otutablelist))
{
  write.csv(x = get(otutablelist[i]), file = paste0(otutablelist[i], "_with_tax.csv"), row.names = FALSE)
}

# write Otutablefull_A_97_minOTU_13 tables to disk

Otutablelist <- ls(pattern = "Otutablefull")
Otutablelist
for(i in 1:length(Otutablelist))
{
  write.csv(x = get(Otutablelist[i]), file = paste0(Otutablelist[i], "_filtered_with_tax.csv"), row.names = FALSE)
}


# save Comm_analysis_list_A1_97_minOTU_13 lists to disk, as RDS objects. Can be read in again via readRDS.

commlist <- ls(pattern = "Comm_analysis") # gets all filenames with "Comm_analysis" in the name
commlist

for(i in 1:length(commlist))
{
	saveRDS(get(commlist[i]), file = paste0(commlist[i], ".rds"))
}

# To read in a list for analyses
community <- readRDS("Comm_analysis_list_A2_97_minOTU_8.rds")
otus <- community[[2]]
env <- community[[1]]




########################################################################
#### Community analyses
########################################################################
# Read in all saved community lists for analysis. This works because the min_OTU threshold was the same for all pools. If not the same, then need to use different filenames

folder <- c("A", "B", "C", "D", "E", "F")

for(i in folder)
{
  for(j in 1:3)
  {
    for(k in 97:97)
    {
      assign(paste0("community", i, j), readRDS(paste0("Comm_analysis_list_",i,j,"_", k, "_minOTU_8.rds")))
    }
  }
}

# do NMDS analysis to see basic patterns ####
bodypart <- setNames(as.data.frame(as.vector(c("body", "leg", "body", "leg", "body", "leg", "body", "leg"))), c("bodypart"))
evenness <- setNames(as.data.frame(as.numeric(c("3", "3", "1", "1", "2", "2", "4", "4"))), c("evenness"))

# folder <- "A"
# pool <- 1

# get(paste0("community","A","3"))[[2]] # this syntax works to extract an object from a list that is constructed from a dynamic variable

# NMDS for all communities
for(i in folder)
{
  for(j in 1:3)
  {
    env <- get(paste0("community", i, j))[[1]]
    env <- bind_cols(env, bodypart, evenness)
    community <- get(paste0("community", i, j))[[2]]
    cat("\n\n", "community", i, j, "\n", sep = "")
    (sprichness <- specnumber(community, groups = env$sample_names, MARGIN = 1)) # number of species per site. NB can't calculate Chao2 because each sample has n=1.
    
    community.jmds <- metaMDS(community, distance = "jaccard", trymax = 40, binary=FALSE)
    # community.jmds <- metaMDS(community, distance = "jaccard", binary = FALSE, previous.best = community.jmds)
    assign(paste0("community", i, j, ".jmds"), community.jmds)
    rm(community.jmds)
    rm(community)
    rm(commlist)
  }
}

# calculate species richness in each sample pool
for(i in folder)
{
  for(j in 1:3)
  {
    env <- get(paste0("community", i, j))[[1]]
    env <- bind_cols(env, bodypart, evenness)
    community <- get(paste0("community", i, j))[[2]]
    cat("sample pool ", i, j, ":  species richness", "\n", sep = "")
    sprichness <- specnumber(community, groups = env$sample_names, MARGIN = 1) # number of species per site. NB can't calculate Chao2 because each sample has n=1.
    cat(sprichness, "\n")
    cat("sample pool ", i, j, ":  Read totals", "\n", sep = "")
    cat(rowSums(community), "\n\n")
    rm(community)
  }
}

# stressplots 
par(mfrow=c(3,3))
for(i in folder[1:3])
{
  for(j in 1:3)
  {
stressplot(get(paste0("community", i, j, ".jmds")), main = paste0("community", i, j, ".jmds"))
  }
}

for(i in folder[4:6])
{
  for(j in 1:3)
  {
    stressplot(get(paste0("community", i, j, ".jmds")), main = paste0("community", i, j, ".jmds"))
  }
}
par(mfrow=c(1,1))

# **compare pools with the same tags at optimal Tm**
# A1:B1, A2:B2, A3:B3
# 
# Have to remove sample mmmmbody (= row 7) from A2 and B2 because of failed PCR of this sample in A2
communityA2_nommmmbody <- communityA2[[2]]
rowSums(communityA2_nommmmbody)
envA2 <- communityA2[[1]]
communityA2_nommmmbody <- communityA2_nommmmbody %>% filter(envA2 != "mmmmbody")
rowSums(communityA2_nommmmbody)
communityB2_nommmmbody <- communityB2[[2]]
envB2 <- communityB2[[1]]
rowSums(communityB2_nommmmbody)
communityB2_nommmmbody <- communityB2_nommmmbody %>% filter(envB2 != "mmmmbody")
rowSums(communityB2_nommmmbody)

# redo NMDS
communityA2_nommmmbody.jmds <- metaMDS(communityA2_nommmmbody, distance = "jaccard", trymax = 40, binary=FALSE)
communityB2_nommmmbody.jmds <- metaMDS(communityB2_nommmmbody, distance = "jaccard", trymax = 20, binary=FALSE)


procrustes(communityA1.jmds, communityB1.jmds)
procrustes(communityA2_nommmmbody.jmds, communityB2_nommmmbody.jmds)
procrustes(communityA3.jmds, communityB3.jmds)
protestA1B1 <- protest(communityA1.jmds, communityB1.jmds)
protestA2B2_nommmmbody <- protest(communityA2_nommmmbody.jmds, communityB2_nommmmbody.jmds)
protestA2B2 <- protest(communityA2.jmds, communityB2.jmds)
protestA3B3 <- protest(communityA3.jmds, communityB3.jmds)

par(mfrow=c(2,2))
plot(protestA1B1, main = "A1 vs B1")
plot(protestA2B2, main = "A2 vs B2")
plot(protestA2B2_nommmmbody, main = "A2 vs B2 (no mmmmbody)")
plot(protestA3B3, main = "A3 vs B3")
par(mfrow=c(1,1))


# **compare pools with different tags at optimal Tm**
# A1n:A2n, A1:A3, A2n:A3n, A1:B2, A1:B3, A2n:B1n, A2n:B3n, A3:B1, A3:B2, B1:B2, B1:B3, B2:B3  # A2n = nommmmbody row

# Remove sample mmmmbody (= row 7) from A1, A3, B1, and B3 because of failed PCR of mmmmbody in A2
# communityA2_nommmmbody already exists

communityA1_nommmmbody <- communityA1[[2]]
rowSums(communityA1_nommmmbody)
envA1 <- communityA1[[1]]
communityA1_nommmmbody <- communityA1_nommmmbody %>% filter(envA1 != "mmmmbody")
rowSums(communityA1_nommmmbody)

communityA3_nommmmbody <- communityA3[[2]]
rowSums(communityA3_nommmmbody)
envA3 <- communityA3[[1]]
communityA3_nommmmbody <- communityA3_nommmmbody %>% filter(envA3 != "mmmmbody")
rowSums(communityA3_nommmmbody)

communityB1_nommmmbody <- communityB1[[2]]
rowSums(communityB1_nommmmbody)
envB1 <- communityB1[[1]]
communityB1_nommmmbody <- communityB1_nommmmbody %>% filter(envB1 != "mmmmbody")
rowSums(communityB1_nommmmbody)

communityB3_nommmmbody <- communityB3[[2]]
rowSums(communityB3_nommmmbody)
envB3 <- communityB3[[1]]
communityB3_nommmmbody <- communityB3_nommmmbody %>% filter(envB3 != "mmmmbody")
rowSums(communityB3_nommmmbody)

# redo NMDS
communityA1_nommmmbody.jmds <- metaMDS(communityA1_nommmmbody, distance = "jaccard", trymax = 40, binary=FALSE)
communityA3_nommmmbody.jmds <- metaMDS(communityA3_nommmmbody, distance = "jaccard", trymax = 20, binary=FALSE)
communityB1_nommmmbody.jmds <- metaMDS(communityB1_nommmmbody, distance = "jaccard", trymax = 20, binary=FALSE)
communityB3_nommmmbody.jmds <- metaMDS(communityB3_nommmmbody, distance = "jaccard", trymax = 20, binary=FALSE)

# A1n:A2n, A1:A3, A2n:A3n, A1:B2, A1:B3, A2n:B1n, A2n:B3n, A3:B1, A3:B2, B1:B2, B1:B3, B2:B3  # A2n = nommmmbody row
procrustes(communityA1_nommmmbody.jmds, communityA2_nommmmbody.jmds)
procrustes(communityA1.jmds, communityA3.jmds)
procrustes(communityA2_nommmmbody.jmds, communityA3_nommmmbody.jmds)
procrustes(communityA1.jmds, communityB2.jmds)
procrustes(communityA1.jmds, communityB3.jmds)
procrustes(communityA2_nommmmbody.jmds, communityB1_nommmmbody.jmds)
procrustes(communityA2_nommmmbody.jmds, communityB3_nommmmbody.jmds)
procrustes(communityA3.jmds, communityB1.jmds)
procrustes(communityA3.jmds, communityB2.jmds)
procrustes(communityB1.jmds, communityB2.jmds)
procrustes(communityB1.jmds, communityB3.jmds)
procrustes(communityB2.jmds, communityB3.jmds)


protestA1A2 <- protest(communityA1_nommmmbody.jmds, communityA2_nommmmbody.jmds)
protestA1A3 <- protest(communityA1.jmds, communityA3.jmds)
protestA2A3 <- protest(communityA2_nommmmbody.jmds, communityA3_nommmmbody.jmds)
protestA1B2 <- protest(communityA1.jmds, communityB2.jmds)
protestA1B3 <- protest(communityA1.jmds, communityB3.jmds)
protestA2B1 <- protest(communityA2_nommmmbody.jmds, communityB1_nommmmbody.jmds)
protestA2B3 <- protest(communityA2_nommmmbody.jmds, communityB3_nommmmbody.jmds)
protestA3B1 <- protest(communityA3.jmds, communityB1.jmds)
protestA3B2 <- protest(communityA3.jmds, communityB2.jmds)
protestB1B2 <- protest(communityB1.jmds, communityB2.jmds)
protestB1B3 <- protest(communityB1.jmds, communityB3.jmds)
protestB2B3 <- protest(communityB2.jmds, communityB3.jmds)


par(mfrow=c(3,4))
# plot(protestA1A2, main = "A1 vs A2")
plot(protestA1A3, main = "A1 vs A3")
# plot(protestA2A3, main = "A2 vs A3")
plot(protestA1B2, main = "A1 vs B2")
plot(protestA1B3, main = "A1 vs B3")
# plot(protestA2B1, main = "A2 vs B1")
# plot(protestA2B3, main = "A2 vs B3")
plot(protestA3B1, main = "A3 vs B1")
plot(protestA3B2, main = "A3 vs B2")
plot(protestB1B2, main = "B1 vs B2")
plot(protestB1B3, main = "B1 vs B3")
plot(protestB2B3, main = "B2 vs B3")
par(mfrow=c(1,1))



# Now do the same as above but with C and D and no removal of mmmmbody


# C1:D1, C2:D2, C3:D3

procrustes(communityC1.jmds, communityD1.jmds)
procrustes(communityC2.jmds, communityD2.jmds)
procrustes(communityC3.jmds, communityD3.jmds)
protestC1D1 <- protest(communityC1.jmds, communityD1.jmds)
protestC2D2 <- protest(communityC2.jmds, communityD2.jmds)
protestC3D3 <- protest(communityC3.jmds, communityD3.jmds)

par(mfrow=c(2,2))
plot(protestC1D1, main = "C1 vs D1")
plot(protestC2D2, main = "C2 vs D2")
plot(protestC3D3, main = "C3 vs D3")
par(mfrow=c(1,1))



# different tags at optimal temperature
# C1:C2, C1:C3, C2:C3, C1:D2, C1:D3, C2:D1, C2:D3, C3:D1, C3:D2, D1:D2, D1:D3, D2:D3

procrustes(communityC1.jmds, communityC2.jmds)
procrustes(communityC1.jmds, communityC3.jmds)
procrustes(communityC2.jmds, communityC3.jmds)
procrustes(communityC1.jmds, communityD2.jmds)
procrustes(communityC1.jmds, communityD3.jmds)
procrustes(communityC2.jmds, communityD1.jmds)
procrustes(communityC2.jmds, communityD3.jmds)
procrustes(communityC3.jmds, communityD1.jmds)
procrustes(communityC3.jmds, communityD2.jmds)
procrustes(communityD1.jmds, communityD2.jmds)
procrustes(communityD1.jmds, communityD3.jmds)
procrustes(communityD2.jmds, communityD3.jmds)

(protestC1C2 <- protest(communityC1.jmds, communityC2.jmds))
(protestC1C3 <- protest(communityC1.jmds, communityC3.jmds))
(protestC2C3 <- protest(communityC2.jmds, communityC3.jmds))
(protestC1D2 <- protest(communityC1.jmds, communityD2.jmds))
(protestC1D3 <- protest(communityC1.jmds, communityD3.jmds))
(protestC2D1 <- protest(communityC2.jmds, communityD1.jmds))
protestC2D3 <- protest(communityC2.jmds, communityD3.jmds)
protestC3D1 <- protest(communityC3.jmds, communityD1.jmds)
protestC3D2 <- protest(communityC3.jmds, communityD2.jmds)
protestD1D2 <- protest(communityD1.jmds, communityD2.jmds)
protestD1D3 <- protest(communityD1.jmds, communityD3.jmds)
protestD2D3 <- protest(communityD2.jmds, communityD3.jmds)


par(mfrow=c(3,4))
plot(protestC1C2, main = "C1 vs C2")
plot(protestC1C3, main = "C1 vs C3")
plot(protestC2C3, main = "C2 vs C3")
plot(protestC1D2, main = "C1 vs D2")
plot(protestC1D3, main = "C1 vs D3")
plot(protestC2D1, main = "C2 vs D1")
plot(protestC2D3, main = "C2 vs D3")
plot(protestC3D1, main = "C3 vs D1")
plot(protestC3D2, main = "C3 vs D2")
plot(protestD1D2, main = "D1 vs D2")
plot(protestD1D3, main = "D1 vs D3")
plot(protestD2D3, main = "D2 vs D3")
par(mfrow=c(1,1))

pairwiseCD <- c("C1C2", "C1C3", "C2C3", "C1D2", "C1D3", "C2D1", "C2D3", "C3D1", "C3D2", "D1D2", "D1D3", "D2D3")

correlationsCD <- vector("double", 12)
for(i in pairwiseCD))
{
  correlationsCD[[i]] <- get(paste0("protest", i))[["scale"]]
}


correlationsCD[[1]] <- protestC1C2[["scale"]]
correlationsCD[[1]]
