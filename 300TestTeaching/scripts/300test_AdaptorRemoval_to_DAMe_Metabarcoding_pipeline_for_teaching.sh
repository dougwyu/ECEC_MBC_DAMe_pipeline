#!/bin/bash
set -e
set -u
set -o pipefail
##########################################################################################################
##########################################################################################################
# script for metabarcoding:  Schirmer et al (2015, NAR) + Zepeda-Mendoza et al (2016) DAMe pipeline
# This version is for teaching, and so only analyses one folder and omits as many bash tricks as possible.
# Most importantly, no loops, to make the code easier to follow.
##########################################################################################################
###########################################################################################################

# Background knowledge
# There is no substitute for first learning about how to use Unix and bash scripting syntax. I found Vince Buffalo's book, Bioinformatic Data Skills to be an excellent introduction. If you are new to Unix, first read Buffalo's Chapters 2, 3, and 7. Then, the downstream analysis is carried out in the R language.  Buffalo's book has an introduction to R, and there are many others online and in books.


####################################################################################################
# Set up folders and pathnames
####################################################################################################
# Set up folders as shown in Readme.md
# Your HOMEFOLDER is the folder that contains the subsidiary folders:  data/, data/seqs/, analysis/, and scripts/
# In this script, the name of the HOMEFOLDER is 300TestTeaching/.
# Thus, run this script from inside ~/PATH_TO/300TestTeaching/scripts/
# The fastq files are placed in 300TestTeaching/data/seqs/
# The DAMe supporting files are placed in 300TestTeaching/data/
# The analysis outputs are placed in in 300TestTeaching/analysis/

# Set path_to_your_homefolder. Do not use the ~ shortcut
HOMEFOLDER="/Users/Negorashi2011/Xiaoyangmiseqdata/MiSeq_20170410/300Test/scripts/300TestTeaching/"
echo "Home folder is" ${HOMEFOLDER}
SEQS="data/seqs/"
ANALYSIS="analysis/"


####################################################################################################
# Install software
####################################################################################################
# The script uses quite a bit of open-source software and is **written for macOS**.  Much of this software can be installed using Homebrew, which is a "package manager" for macOS.  Install that, and most software can then be installed and managed with a single command.
# In principle, you only need to do this once

## Atom.io:  Text editor
# Download from https://atom.io
# Install platformio-ide-terminal, which allows you to send commands to the terminal within Atom
# In Atom, navigate to Atom/Preferences.../+Install and type platformio-ide-terminal into the "Search packages" window and press Enter. Click on the blue Install button.
# A new terminal can be opened in the bottom left of the Atom window, by clicking on the + sign.

## Homebrew
# go to http://brew.sh and follow the instructions for installing Homebrew on macOS
# after Homebrew is installed, run these brew installations
brew tap homebrew/science # a "tap" is a source of "installation formulae" of specialist software, here, bioinformatics
brew install git
brew install coreutils
brew install gnu-sed # (gsed == GNU version of sed == Linux version of sed)
brew install adapterremoval # https://github.com/MikkelSchubert/adapterremoval
brew install sickle # https://github.com/najoshi/sickle
brew install homebrew/science/spades # http://cab.spbu.ru/software/spades/
brew install homebrew/science/pandaseq # https://github.com/neufeld/pandaseq/releases
brew install homebrew/science/vsearch # https://github.com/torognes/vsearch
brew install seqtk # https://github.com/lh3/seqtk

## DAMe
cd /usr/local/bin
git clone https://github.com/shyamsg/DAMe.git # can read tags with heterogeneity spacers
# no building needed

## RStudio
# install from https://www.rstudio.com

## R
# install from https://cran.r-project.org
# Then install R packages by running these commands:
install.packages(c("tidyverse", "data.table", "vegan", "car", "RColorBrewer"), dependencies = TRUE)
source("https://bioconductor.org/biocLite.R") # to install bioinformatics packages
biocLite("phyloseq") # install phyloseq

## Sumatra
cd ~/Desktop/
wget https://git.metabarcoding.org/obitools/sumatra/uploads/251020bbbd6c6595cb9fce6077e29952/sumatra_v1.0.20.tar.gz
tar -zxvf sumatra_v1.0.20.tar.gz
cd sumatra_v1.0.20/
make CC=clang # disables OpenMP, which isn't on macOS
mv sumatra /usr/local/bin

## Sumaclust
cd ~/Desktop/
wget https://git.metabarcoding.org/obitools/sumaclust/uploads/69f757c42f2cd45212c587e87c75a00f/sumaclust_v1.0.20.tar.gz
tar -zxvf sumaclust_v1.0.20.tar.gz
cd sumaclust_v1.0.20/
make CC=clang # disables OpenMP, which isn't on macOS
mv sumaclust /usr/local/bin/


####################################################################################################
# Set variable values
####################################################################################################
DAME="/usr/local/bin/DAMe/bin/" # this should be the pathname for the binaries on DAMe
echo "The DAMe binaries are in" ${DAME}
SUMASIM=97 # sumaclust similarity in percentage (0-100, e.g. 97)
echo "Sumaclust similarity percentage is" ${SUMASIM}"%."
MINPCR=2 # min contig length after assembly if i'm not running as a bash script
MINREADS=4 # min contig length after assembly if i'm not running as a bash script
MINPCR_1=1 # for RSI analysis and to run without filtering by min PCR number or min copy number per pcr
MINREADS_1=1 # for RSI analysis and to run without filtering by min PCR number or min copy number per pcr
MINLEN=300 # minimum length of a read
MAXLEN=320 # maximum length of a read
POOLS=3 # number of times that a library was PCRd (following DAMe protocol, using different twin-tags for each PCR rxn)
ARTHMINPROB=0.8 # lowest allowable RDP-Classifier assignment probability


# For teaching, I reduced the original fastq files down to 10%.


####################################################################################################
####################################################################################################
# STAGE ONE:  Trimming, Denoising, and Merging the raw sequencing files
####################################################################################################
####################################################################################################

cd ${HOMEFOLDER}${SEQS} # cd into the seqs/ folder to start stage one.
touch software_versions.txt # to store software versions


####################################################################################################
# 1.1 Use AdapterRemoval to trim Illumina sequencing adapters. (there are a few adapters still left in the raw data)
####################################################################################################

AdapterRemoval -h
AdapterRemoval -v  # confirm the version
echo "AdapterRemoval ver. 2.2.2" >> software_versions.txt  # cannot redirect the output of AdapterRemoval -v

# identifies whether and which Illumina adapters are left in the raw data. Only needs to be run for one library (e.g. B1)
AdapterRemoval --identify-adapters --file1 B1_S4_L001_R1_001_0.1.fastq.gz --file2 B1_S4_L001_R2_001_0.1.fastq.gz

# Teaching note. IMPORTANT:  Before running the following commands, note that if the command wraps around to a 2nd line, pressing ctrl-enter will move the cursor into the second line of the current command, NOT into the next command. If you do not move your cursor into the next command (using the down-arrow key), pressing ctrl-enter again will merely send the same command again. Usually, this just wastes time, but it could cause serious errors, if the command is not written properly. Always check where your cursor is before pressing ctrl-enter.

# removes adapters from the fastq reads and outputs new files with fixed fastq files.
AdapterRemoval  --file1 B1_S4_L001_R1_001_0.1.fastq.gz --file2 B1_S4_L001_R2_001_0.1.fastq.gz --basename AdaptermvB1 --threads 3 --trimns
AdapterRemoval  --file1 B2_S5_L001_R1_001_0.1.fastq.gz --file2 B2_S5_L001_R2_001_0.1.fastq.gz --basename AdaptermvB2 --threads 3 --trimns
AdapterRemoval  --file1 B3_S6_L001_R1_001_0.1.fastq.gz --file2 B3_S6_L001_R2_001_0.1.fastq.gz --basename AdaptermvB3 --threads 3 --trimns

# Teaching note. We are applying the same command (AdapterRemoval) to multiple files. Normally, instead of writing the command for each set of files, we would write a loop that inserts the correct filename. This is a less buggy approach but is harder to understand in the beginning. The main risk is that I use the wrong filename (e.g. B1 instead of B2 in the second command. This kind of editing mistake is easy to make.)


####################################################################################################
# 1.2 Use Sickle to trim away low-quality nucleotides from the 3' ends of the reads (does not try to fix errors)
####################################################################################################

sickle -h
sickle --version >> software_versions.txt # sickle version 1.33

# Usage: sickle pe [options] -f <paired-end forward fastq file> -r <paired-end reverse fastq file> -o <trimmed PE forward file> -p <trimmed PE reverse file> -s <trimmed singles file> -t <quality type>

# The sickle_Single_B1.fq files are single (no longer paired) reads. These will not be used in the DAMe pipeline, as we are using twin tags.
sickle pe -f AdaptermvB1.pair1.truncated -r AdaptermvB1.pair2.truncated -o sickle_B1_R1.fq -p sickle_B1_R2.fq -s sickle_Single_B1.fq -t sanger > sickleB1.out
sickle pe -f AdaptermvB2.pair1.truncated -r AdaptermvB2.pair2.truncated -o sickle_B2_R1.fq -p sickle_B2_R2.fq -s sickle_Single_B2.fq -t sanger > sickleB2.out
sickle pe -f AdaptermvB3.pair1.truncated -r AdaptermvB3.pair2.truncated -o sickle_B3_R1.fq -p sickle_B3_R2.fq -s sickle_Single_B3.fq -t sanger > sickleB3.out

# Now you can remove the AdapterRemoval output files
rm Adaptermv*.*


####################################################################################################
# 1.3 Use SPAdes to correct sequencing and PCR errors, via BayesHammer
####################################################################################################

spades.py -h
spades.py -v # confirm the version
echo "SPAdes v3.10.1" >> software_versions.txt # SPAdes v3.10.1

# Each pair of (the reduced) fastq files typically requires 4-5 minutes. The original files required ~45 minutes per library, on an Core i5 MacBook Pro, using only 3 threads. We are not denoising the unpaired reads.
spades.py --only-error-correction -1 sickle_B1_R1.fq -2 sickle_B1_R2.fq -o SPAdes_hammer_B1 -t 3
spades.py --only-error-correction -1 sickle_B2_R1.fq -2 sickle_B2_R2.fq -o SPAdes_hammer_B2 -t 3
spades.py --only-error-correction -1 sickle_B3_R1.fq -2 sickle_B3_R2.fq -o SPAdes_hammer_B3 -t 3

# Now you can remove the sickle files
rm sickle_*.fq


####################################################################################################
# 1.4 Use PandaSeq to merge the paired reads
####################################################################################################

pandaseq -h
pandaseq -v # pandaseq 2.11
echo "pandaseq 2.11" >> software_versions.txt

# Usage: pandaseq -f forward.fastq -i index.fastq -r reverse.fastq [-6] [-A algorithm:parameters] [-B] [-C filter] [-D threshold] [-F] [-G log.txt.bz2] [-L length] [-N] [-O length] [-T threads] [-U unaligned.txt] [-W output.fasta.bz2] [-a] [-d flags] [-g log.txt] [-h] [-j] [-k kmers] [-l length] [-o length] [-p primer] [-q primer] [-t threshold] [-u unaligned.txt] [-v] [-w output.fasta]

# Each pair of fastq files takes ~ 2 secs
pandaseq -f SPAdes_hammer_B1/corrected/sickle_B1_R1.00.0_0.cor.fastq.gz -r SPAdes_hammer_B1/corrected/sickle_B1_R2.00.0_0.cor.fastq.gz -A simple_bayesian -d bfsrk -F -N -T 7 -g pandaseq_log_B1.txt -w sickle_cor_panda_B1.fastq

pandaseq -f SPAdes_hammer_B2/corrected/sickle_B2_R1.00.0_0.cor.fastq.gz -r SPAdes_hammer_B2/corrected/sickle_B2_R2.00.0_0.cor.fastq.gz -A simple_bayesian -d bfsrk -F -N -T 7 -g pandaseq_log_B2.txt -w sickle_cor_panda_B2.fastq

pandaseq -f SPAdes_hammer_B3/corrected/sickle_B3_R1.00.0_0.cor.fastq.gz -r SPAdes_hammer_B3/corrected/sickle_B3_R2.00.0_0.cor.fastq.gz -A simple_bayesian -d bfsrk -F -N -T 7 -g pandaseq_log_B3.txt -w sickle_cor_panda_B3.fastq


####################################################################################################
# 1.5 gzip the final fastq files and remove working files
####################################################################################################

echo "B1"
gzip sickle_cor_panda_B1.fastq
echo "B2"
gzip sickle_cor_panda_B2.fastq
echo "B3"
gzip sickle_cor_panda_B3.fastq

# Now you can remove the SPAdes output
rm -rf SPAdes_hammer*/
# remove pandaseq_log_txt files
rm pandaseq_log_*.txt

# The sickle_cor_panda_B{1,2,3} are your final fastq files.
# They have been adapter-trimmed, quality-trimmed, denoised, and merged.
# Three pairs of fastq files (forward and reverse reads) have been reduced to three fastq files (merged reads).
# These are the input files for DAMe.
# You will now use DAMe to look for reads that appear in more than one library (B1, B2, B3).



####################################################################################################
####################################################################################################
# STAGE TWO:  DAMe filtering the merged reads to keep only the sequences that appear in
#             multiple PCRs and in multiple copies
####################################################################################################
####################################################################################################


####################################################################################################
# 2.1 Place the libraries in a folder for that experiment (B) and then into separate folders
####################################################################################################

mkdir folder_B/ # make a folder for the experiment (here, B)
mkdir folder_B/pool1/
mkdir folder_B/pool2/
mkdir folder_B/pool3/

mv sickle_cor_panda_B1.fastq.gz folder_B/pool1/
mv sickle_cor_panda_B2.fastq.gz folder_B/pool2/
mv sickle_cor_panda_B3.fastq.gz folder_B/pool3/


####################################################################################################
# 2.2 Sort through each fastq file and determine how many of each tag pair is in each fastq file. Each fastq file takes < 1 min
####################################################################################################

${DAME}sort.py -h #
# DAMe does not have version numbers. We are using the branch from https://github.com/shyamsg/DAMe.git
echo "https://github.com/shyamsg/DAMe.git Downloaded on 14 July 2017" >> ${HOMEFOLDER}data/seqs/software_versions.txt

cd ${HOMEFOLDER}data/seqs/folder_B/pool1/
python ${DAME}sort.py -fq sickle_cor_panda_B1.fastq.gz -p ${HOMEFOLDER}data/Primers_COILeray.txt -t ${HOMEFOLDER}data/Tags_300test_COI.txt

cd ${HOMEFOLDER}data/seqs/folder_B/pool2/
python ${DAME}sort.py -fq sickle_cor_panda_B2.fastq.gz -p ${HOMEFOLDER}data/Primers_COILeray.txt -t ${HOMEFOLDER}data/Tags_300test_COI.txt

cd ${HOMEFOLDER}data/seqs/folder_B/pool3/
python ${DAME}sort.py -fq sickle_cor_panda_B3.fastq.gz -p ${HOMEFOLDER}data/Primers_COILeray.txt -t ${HOMEFOLDER}data/Tags_300test_COI.txt

# Teaching note.  Look at the sort.py output for B3:

     # Number of erroneous sequences in file sickle_cor_panda_B3.fastq.gz (with errors in the sequence of primer or tags, or no barcode amplified):
     # 12033
     #
     #     No sequence between primers         : 0
     #     Tags pair not found                 : 2024
     #     F primer found, R' primer not found : 4186
     #     R primer found, F' primer not found : 1268
     #     Neither F nor R primer found        : 4555
     #
     # Number of valid tag pairs found         : 33468
     #     F-R' barcodes found                 : 25153
     #     R-F' barcodes found                 : 8315


# sort.py looks for the tag pairs that you expect to be in this experiment B.
# The expected tags are listed in PSinfo_300test_COIB.txt. Note that sort.py looks for all tags from
# all 3 PCRs.
# There are 12033 reads that have tag errors (e.g. missing primer or tags) and 33468 sequences with valid tag pairs.
# This version of DAMe does not accept primers or tags with any mismatches, so we are rejecting reads that might be useful. There is a plan to add the ability to accept primers and/or tags with n mismatches.


####################################################################################################
# 2.3 Make an overview of the tag combinations.
####################################################################################################

# In the output file, the term "was used" means that the tag pair is in the PSInfo file (i.e. The tag pair was used in PCR).

python ${DAME}splitSummaryByPSInfo.py -h

cd ${HOMEFOLDER}data/seqs/folder_B

python ${DAME}splitSummaryByPSInfo.py -p ${HOMEFOLDER}data/PSinfo_300test_COIB.txt -l 1 -s pool1/SummaryCounts.txt -o splitSummaryByPSInfo_B1.txt
python ${DAME}splitSummaryByPSInfo.py -p ${HOMEFOLDER}data/PSinfo_300test_COIB.txt -l 2 -s pool2/SummaryCounts.txt -o splitSummaryByPSInfo_B2.txt
python ${DAME}splitSummaryByPSInfo.py -p ${HOMEFOLDER}data/PSinfo_300test_COIB.txt -l 3 -s pool3/SummaryCounts.txt -o splitSummaryByPSInfo_B3.txt

# Teaching Note.  Look at the syntax of splitSummaryByPSInfo.py.  -l 1 indicates that this is for "pool 1", which is the first PCR replicate of experiment B.  DAMe uses this information to search only for tag pairs used in pool 1, in the file PSinfo_300test_COIB.txt.

# Look at the output file splitSummaryByPSInfo_B1.txt.  In the first section ("Tag combinations where the tag pair was used"), Tag13-Tag13 to Tag4-Tag4 and also Tag7-Tag7 have many reads assigned. This is good because these tag pairs were used to amplify insect soups. In contrast, Tag10-Tag10, Tag11-Tag11, Tag12-Tag12 are the extraction blanks. There are very few reads in them. Tag8-Tag8 is the pooled PCR blank and also has few reads.

# Look at the output file splitSummaryByPSInfo_B3.txt.  There are fewer rows in the first section. The missing tag pairs have zero reads assigned to them, and they were extraction blanks. So this PCR reaction was very clean.


####################################################################################################
# 2.4 Generate a heatmap of the tag pair read information.
####################################################################################################

# This step requires R (and helpful to have RStudio)
# Installing R should also have installed the binary "Rscript", which we will use to invoke R from the command line

cd ${HOMEFOLDER}data/seqs/folder_B/pool1/
Rscript --vanilla --verbose ${HOMEFOLDER}scripts/heatmap_for_teaching.R

cd ${HOMEFOLDER}data/seqs/folder_B/pool2/
Rscript --vanilla --verbose ${HOMEFOLDER}scripts/heatmap_for_teaching.R

cd ${HOMEFOLDER}data/seqs/folder_B/pool3/
Rscript --vanilla --verbose ${HOMEFOLDER}scripts/heatmap_for_teaching.R

# then move heatmaps from inside the pool folders into folder_B/
cd ${HOMEFOLDER}data/seqs/folder_B/
mv pool1/heatmap.pdf heatmap_B1.pdf
mv pool2/heatmap.pdf heatmap_B2.pdf
mv pool3/heatmap.pdf heatmap_B3.pdf

# Teaching note.  Examine the heatmaps. The pink squares should be on the diagonal only. The blue squares indicate low numbers of reads of non-twin tag pairs. These are tag jumps.

# Sometimes, a PCR-setup error can use the wrong (e.g. non-twin) tag pair for a sample. The heatmaps can find and rescue these samples. See the DAMe paper for an example.


####################################################################################################
# 2.5 Check if the PCRs were carried out correctly.
####################################################################################################

# Because we independently carried out PCR 3 times for each sample, the 3 PCRs for each sample should show high similarity. If not, then something went wrong in the PCR setup (e.g. PCR failure, use of wrong sample). This step calculates the Renkonen Similarity Index (RSI).  See the DAMe paper for explanation.

# First step is to use filter.py but keep *all* sequences and then run RSI.py. This step can take > 1 hr for a typical sized file (3 pools in )

# The folders holding the three pools in Experiment B need to be named pool1/, pool2/, and pool3/.
# Run filter.py inside folder_B/

cd ${HOMEFOLDER}data/seqs/folder_B/

mkdir Filter_min${MINPCR_1}PCRs_min${MINREADS_1}copies_B # make a folder for the output

python ${DAME}filter.py -psInfo ${HOMEFOLDER}data/PSinfo_300test_COIB.txt -x ${POOLS} -y ${MINPCR_1} -p ${POOLS} -t ${MINREADS_1} -l ${MINLEN} -o Filter_min${MINPCR_1}PCRs_min${MINREADS_1}copies_B

python ${DAME}RSI.py --explicit -o RSI_output_B_min${MINPCR_1}PCR_min${MINREADS_1}copy.txt Filter_min${MINPCR_1}PCRs_min${MINREADS_1}copies_B/Comparisons_${POOLS}PCRs.txt

# Teaching note.  The RSI.py script compares the similarity of all the PCRs (pools) within each sample. For example, sample mmmmbody was PCRd three times, and the three PCRs are only moderately dissimilar from each other (0.4 < RSI < 0.6), which is acceptable.  Note that the xb samples are completely dissimilar (RSI = 1 = have zero sequences in common). This is expected for these very small negative controls and is a good result.  For more detail, see the DAMe paper.

# Sample	ReplicateA	ReplicateB	RSI
# mmmmbody	1	2	0.527783258807
# mmmmbody	1	3	0.574232901387
# mmmmbody	2	3	0.560717849504
# Hhmlbody	1	2	0.532822968348
# Hhmlbody	1	3	0.537849521358
# Hhmlbody	2	3	0.52770501587
# hlllleg	1	2	0.525170721128
# hlllleg	1	3	0.537077518076
# hlllleg	2	3	0.499159678843
# xb1	1	2	1.0
# xb1	1	3	1.0
# xb1	2	3	1.0
# xb3	1	2	1.0
# xb3	1	3	1.0
# xb3	2	3	1.0


####################################################################################################
# 2.6 Choose thresholds for filtering by minimum number of PCRS (pools) and copies.
####################################################################################################

python ${DAME}plotLengthFreqMetrics_perSample.py -h

cd ${HOMEFOLDER}data/seqs/folder_B/Filter_min${MINPCR_1}PCRs_min${MINREADS_1}copies_B/

# -n 3 is the number of pools
python ${DAME}plotLengthFreqMetrics_perSample.py -f FilteredReads.fna -p ${HOMEFOLDER}data/PSinfo_300test_COIB.txt -n 3

# Teaching note.  I see that there are a few reads > 320 nucleotides. These will need to be filtered out later if filter.py does not do it automatically in the next step.


####################################################################################################
# 2.7 Filter the reads by presence in a user-defined minimum number of PCRS (pools) and copies.
####################################################################################################

python ${DAME}filter.py -h

# After consideration of the negative and positive controls, splitSummaryByPSInfo, and heatmaps, choose thresholds for filtering out reads.

# The two thresholds are the minimum number of PCRs (out of 3) in which a read must appear (minPCR) and the minimum number of copies of that read in each PCR (minReads). Note that minReads is per PCR, so if minPCR=2 and minReads=4, then the each read must appear at least 2 X 4 = 8 times, at least 4 times in each PCR, not 4 times total. Appearing 8 in one PCR and 0 in another is not acceptable in this scenario.

# You will probably want to run a few times through the OTU clustering step, using different values of minPCR and minReads to find values that result in the lowest number of OTUs from the extraction blank and PCR blank negative controls and highest recovery of OTUs from your positive controls.

# This step runs quite slowly.  ~1.3 hrs per library on my MacBook Pro, but filter.py runs single-threaded, so you can launch multiple threads and run in parallel.  I use GNU parallel to achieve this.

#### If I want to re-run filter.py with different thresholds, I can set new values of MINPCR & MINREADS here.
# MINPCR=2 # min number of PCRs that a sequence has to appear in
# MINREADS=4 # min number of copies per sequence per PCR
# confirm MINPCR and MINREADS values
echo "Each (unique) sequence must appear in at least ${MINPCR} PCRs, with at least ${MINREADS} reads per PCR."

cd ${HOMEFOLDER}data/seqs/folder_B/
# make folder for output
mkdir Filter_min${MINPCR}PCRs_min${MINREADS}copies_B/

python ${DAME}filter.py -psInfo ${HOMEFOLDER}data/PSinfo_300test_COIB.txt -x ${POOLS} -y ${MINPCR} -p ${POOLS} -t ${MINREADS} -l ${MINLEN} -o Filter_min${MINPCR}PCRs_min${MINREADS}copies_B

# Teaching note. The output file is called FilteredReads.fna, inside the Filter_min2PCRs_min4copies_B folder.  This is the product of all your hard work!  The next step is to figure out how best to cluster these reads into OTUs.

# re-run plotLengthFreqMetrics_perSample.py
cd ${HOMEFOLDER}data/seqs/folder_B/Filter_min${MINPCR}PCRs_min${MINREADS}copies_B
python ${DAME}plotLengthFreqMetrics_perSample.py -f FilteredReads.fna -p ${HOMEFOLDER}data/PSinfo_300test_COIB.txt -n 3

# The >320 bp reads are gone.


####################################################################################################
# 2.8 Prepare FilteredReads.fna file for OTU clustering and remove chimeras using vsearch
####################################################################################################

python ${DAME}convertToUSearch.py -h

# MINPCR=2 # these commands are to make it possible to process multiple filter.py outputs, which are saved in different Filter_min folders
# MINREADS=4

# vsearch requires a few seconds per library when applied to FilteredReads.forsumaclust.fna
cd ${HOMEFOLDER}data/seqs/folder_B/Filter_min${MINPCR}PCRs_min${MINREADS}copies_B

python ${DAME}convertToUSearch.py -i FilteredReads.fna -lmin ${MINLEN} -lmax ${MAXLEN}
gsed 's/ count/;size/' FilteredReads.forsumaclust.fna > FilteredReads.forvsearch.fna

vsearch --sortbysize FilteredReads.forvsearch.fna --output FilteredReads.forvsearch_sorted.fna
vsearch --uchime_denovo FilteredReads.forvsearch_sorted.fna --nonchimeras FilteredReads.forvsearch_sorted_nochimeras.fna

gsed 's/;size/ count/' FilteredReads.forvsearch_sorted_nochimeras.fna > FilteredReads.forsumaclust.nochimeras.fna

# remove vsearch uchime working files
rm FilteredReads.forsumaclust.fna
rm FilteredReads.forvsearch.fna
rm FilteredReads.forvsearch_sorted.fna
rm FilteredReads.forvsearch_sorted_nochimeras.fna

# After DAMe filtering, vsearch -uchime_denovo will find <0.1% of reads to be chimeras, which doesn't seem worth doing, but each of these reads could turn into a separate OTU, adding dozens or more extra (but small) OTUs to the final dataset.
# We do not detect chimeras in this dataset because it is so small and because DAMe filtering already does a good job of removing artefact reads.


####################################################################################################
# 2.8 Analyse FilteredReads.fna for pairwise similarities to choose similarity threshold for Sumaclust
####################################################################################################

python ${DAME}assessClusteringParameters.py -h

# Assessing clustering parameters. This takes quite a long time on large fna files and is only feasible on FilteredReads.fna that have been filtered to a few MB.  Takes about 5 mins on my MacBook Pro
# MINPCR=2 # these commands are to make it possible to process multiple filter.py outputs, which are saved in different Filter_min folders
# MINREADS=4
echo "Analysing filter.py output where each (unique) sequence appeared in ≥ ${MINPCR} PCRs, with ≥ ${MINREADS} reads per PCR."

cd ${HOMEFOLDER}data/seqs/folder_B/Filter_min${MINPCR}PCRs_min${MINREADS}copies_B

python ${DAME}assessClusteringParameters.py -i FilteredReads.forsumaclust.nochimeras.fna -mint 0.8 -minR 0.6 -step 0.05 -t 4 -o COIclusterassess_mint08_minR06_step005_Filter_min${MINPCR}PCRs_min${MINREADS}copies_${sample}.pdf

# Teaching note.  The key result is that very few to no reads are similar at 96-97%.  This similarity threshold represents the famous 'barcoding gap,' and we choose 97% as the clustering threshold for Sumaclust

####################################################################################################
# 2.9 Sumaclust clustering and convert Sumaclust output to table format
####################################################################################################

python ${DAME}tabulateSumaclust.py -h

# 97% sumaclust
echo ${SUMASIM} # confirm that there is a similarity value chosen
SUMASIM=97 # if there is no SUMASIM value
echo ${SUMASIM} # confirm that there is a similarity value chosen
cd ${HOMEFOLDER}data/seqs/folder${sample}/Filter_min${MINPCR}PCRs_min${MINREADS}copies_${sample}

# cluster reads into OTUs. Output file assigns each read to a cluster.
sumaclust -t .${SUMASIM} -e FilteredReads.forsumaclust.nochimeras.fna > OTUs_${SUMASIM}_sumaclust.fna

# create OTU table and create a fasta file for OTU representative sequences
python ${DAME}tabulateSumaclust.py -i OTUs_${SUMASIM}_sumaclust.fna -o table_300test_B_${SUMASIM}.txt -blast

mv table_300test_B_${SUMASIM}.txt.blast.txt table_300test_B_${SUMASIM}.fas # chg filename suffix to *.fas


####################################################################################################
# 2.10 Move Sumaclust outputs to ${HOMEFOLDER}/analysis
####################################################################################################

# We now have our OTUs table and OTU representative sequences.
# For downstream analysis, I suggest moving these outputs to the analysis folder

cd ${HOMEFOLDER}${ANALYSIS}
# make folders to hold results.
mkdir OTU_transient_results/
mkdir OTU_transient_results/OTU_tables/

mv ${HOMEFOLDER}${SEQS}folder_B/Filter_min${MINPCR}PCRs_min${MINREADS}copies_B ${HOMEFOLDER}${ANALYSIS}OTU_transient_results/Filter_min${MINPCR}PCRs_min${MINREADS}copies_B
mv ${HOMEFOLDER}${SEQS}software_versions.txt ${HOMEFOLDER}${ANALYSIS}OTU_transient_results/

cd ${HOMEFOLDER}${ANALYSIS}OTU_transient_results/Filter_min${MINPCR}PCRs_min${MINREADS}copies_B

rm Comparisons_${POOLS}PCRs.fasta # large file that doesn't need to be kept
rm Comparisons_${POOLS}PCRs.txt # large file that doesn't need to be kept

# I find it useful to make a copy of the OTU table and representative sequences into a separate folder.
# I work with the files in this folder.
cp table_300test_B_*.txt ${HOMEFOLDER}${ANALYSIS}OTU_transient_results/OTU_tables/ # copy OTU tables to OTU_tables folder
cp table_300test_B_*.fas ${HOMEFOLDER}${ANALYSIS}OTU_transient_results/OTU_tables/ # copy OTU repr seqs to OTU_tables folder

# change name of OTU_transient_results folder to include filter.py thresholds and timestamp
mv ${HOMEFOLDER}${ANALYSIS}OTU_transient_results/ ${HOMEFOLDER}${ANALYSIS}OTUs_min${MINPCR}PCRs_min${MINREADS}copies_"$(date +%F_time-%H%M)"/


####################################################################################################
# 2.11 Assign taxonomies to the OTU representative sequences, using RDP Classifier and the Midori dataset
####################################################################################################

# I upload the OTU fasta files to a supercluster (hpc.uea.ac.uk) and assign taxonomies via RDP Classifier on the Midori database. I had previously trained Classifier on Midori. I haven't tried running on my macOS because on the cluster, this step requires ~18GB of RAM to run and crashes on the cluster if it doesn't get that.  RDP classification takes ~4 mins with this dataset.

# bsub file
####################################################################################################
#!/bin/sh
#BSUB -q mellanox-ib     # on hpc.uea.ac.uk, mellanox-ib (168 hours = 7 days)
#BSUB -J RDPclmidAF
#BSUB -oo RDPclmidAF.out
#BSUB -eo RDPclmidAF.err
#BSUB -R "rusage[mem=64000]"  # max mem for mellanox-ib is 128GB = 128000
#BSUB -M 64000
#BSUB -B        # sends email to me when job starts
#BSUB -N        # sends email to me when job finishes
# . /etc/profile
# module load java/jdk1.8.0_51
# java -Xmx64g -jar ~/scripts/RDPTools/classifier.jar classify -t ~/midori/RDP_out/rRNAClassifier.properties -o ~/Yahan/300Test/teaching/table_300test_B_97.RDPmidori.txt ~/Yahan/300Test/teaching/table_300test_B_97.fas
####################################################################################################

# Download the RDP output file (table_300test_B_97.RDPmidori.txt) back to the OTU_tables folder and filter out the OTUs that were not assigned to Arthropoda with probability >= 0.80 (ARTHMINPROB).  I had previoulsy determined using a leave-one-out test that assignments to phylum (Arthropoda) and class (Insecta, Arachnida) ranks are reliable if the probability >= 0.80.

echo "Assignment probability minimum set to: " ${ARTHMINPROB}

OTUTABLEFOLDER="OTUs_min2PCRs_min4copies_2017-08-17_time-1634"

# We have to use some bash tricks to do this efficiently.   If gsed does not run, try replacing it with sed.

cd ${HOMEFOLDER}${ANALYSIS}/${OTUTABLEFOLDER}/OTU_tables

awk -v arthmin=${ARTHMINPROB} '$8 ~ /Arthropoda/ && $10 >= arthmin { print }' table_300test_B_${SUMASIM}.RDPmidori.txt > table_300test_B_${SUMASIM}.RDPmidori_Arthropoda.txt

gsed -E 's/\t\t/\t/' table_300test_B_${SUMASIM}.RDPmidori_Arthropoda.txt > table_300test_B_${SUMASIM}.RDPmidori_Arthropoda_nodbltab.txt

mv table_300test_B_${SUMASIM}.RDPmidori_Arthropoda_nodbltab.txt table_300test_B_${SUMASIM}.RDPmidori_Arthropoda.txt

seqtk subseq table_300test_B_${SUMASIM}.fas <(cut -f 1 table_300test_B_${SUMASIM}.RDPmidori_Arthropoda.txt) > table_300test_B_${SUMASIM}_Arthropoda.fas
     # <(cut -f 1 table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt) produces the list of sequences to keep

echo "OTU tables"
wc -l table_300test_B_${SUMASIM}.RDPmidori.txt
wc -l table_300test_B_${SUMASIM}.RDPmidori_Arthropoda.txt
echo "fasta files"
grep ">" table_300test_B_${SUMASIM}.fas | wc -l
grep ">" table_300test_B_${SUMASIM}_Arthropoda.fas | wc -l

# The number of rows in the final OTU table should be equal to the number of sequences in the final fasta file:  220.

# 220 OTU sequences is the final product. These OTUs have been filtered for reliability and taxonomic assignment to Arthropoda.

# The third stage is to read the files into R to perform the final filtering steps. In R, we filter out "small" OTUs using the phyloseq method and we filter out contaminant sequences via direct observation of the negative controls and positive controls. We also remove the positive control sequences from the dataset.

# The R script is:  300Test_filtered_OTU_table_for_teaching.R

# The final outputs are a fasta file with the OTU representative sequences, an OTU table (OTU x sample) with taxonomic assignments, and community analysis files (sample X OTU + an environment table (sample X environmental variables)). We can carry out community analysis with the last output.

#### End script here.
exit


####################################################################################################
####################################################################################################
# APPENDIX.  Fixing the RDP Classifier taxonomic assignments by adding in missing ranks
####################################################################################################
####################################################################################################

# NCBI records sometimes omit some of the ranks.  E.g. A sequence might have an identification to species but not have the family. This causes the table of taxonomic assignments to lack sufficient columns, and that prevents R from importing. This is example code to fix that problem on a record by record basis.  Usually, there are few of these problems.

# Hemiptera;Caliscelidae;Bruchomorpha is missing its family in the Midori database. So the identifcation omits the family name. This prevents R from inputting the OTU table.

# cd ${HOMEFOLDER}${ANALYSIS}/${OTUTABLEFOLDER}/OTU_tables
#
# gsed -E 's/Bruchomorpha\tgenus/Caliscelidae\tfamily\t0.50\t\Bruchomorpha\tgenus/' table_300test_B_${SUMASIM}.RDPmidori_Arthropoda.txt > table_300test_B_${SUMASIM}.RDPmidori_Arthropoda_family.txt
#
# mv table_300test_B_${SUMASIM}.RDPmidori_Arthropoda_family.txt table_300test_B_${SUMASIM}.RDPmidori_Arthropoda.txt
#
# grep "Bruchomorpha" table_300test_B_${SUMASIM}.RDPmidori_Arthropoda.txt
#


# DANGEROUS CODE:  RUN ONLY ONCE AFTER GENERATING THE RDP ARTHROPODA-ONLY TABLES, BECAUSE IF RUN MORE THAN ONCE, WILL INSERT THE NEW TAXONOMIC RANK (e.g. Caliscelidae family 0.5)  MORE THAN ONCE
# OTUTABLEFOLDER="OTUs_min2PCRs_min4copies_2017-08-09_time-1249"
# for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
# do
#      for sim in `seq 96 97`  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
#      do
#      done
# done
