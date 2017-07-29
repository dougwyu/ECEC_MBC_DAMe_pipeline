#!/bin/bash
set -e
set -u
set -o pipefail
#######################################################################################
#######################################################################################
# a shell script for metabarcoding:  Schirmer + DAMe pipeline
#######################################################################################
#######################################################################################

# To do
              # add a chimera checking step after I get the OTU tables:  the OTU table gives me the OTU size information (sum of each row's read numbers), and I can use the OTU representative fasta file as the input. This needs a bit of programming. Must use usearch9 uchime2_denovo, because usearch10 uchime3_denovo isn't working. After running, need a step to edit the OTU table and OTU fasta file

              # try bfc instead of spade.py --error_correction-only
              # understand output of assessClusteringParameters.py

              ## some additional DAMe commands to consider (unfold to see details)


# Usage: bash 300test_AdaptorRMV+Panda+DAMe.sh SUMASIM
# e.g. 300test_AdaptorRMV+Panda+DAMe.sh 2 3 97 # OTU must appear in at least 2 PCRs and 3 reads per PCR, and will be clustered at 97%

PIPESTART=$(date)

# run script from ~/PATH_TO/300Test/scripts
# fastq files in 300Test/data/seqs/
# DAMe files in 300Test/data/
# analysis outputs in 300Test/analysis/
# this script in 300Test/scripts/

# set variables
# SUMASIM=$1 # sumaclust similarity in percentage (0-100, e.g. 97)
SUMASIM=97 # sumaclust similarity in percentage (0-100, e.g. 97)
MINPCR=2 # min contig length after assembly if i'm not running as a bash script
MINREADS=3 # min contig length after assembly if i'm not running as a bash script
MINLEN=300
PCRRXNS=3
POOLS=3
echo "Sumaclust similarity percentage is "$SUMASIM"%."
HOMEFOLDER="/Users/Negorashi2011/Xiaoyangmiseqdata/MiSeq_20170410/300Test/"  # do not have a ~ in the path
echo "Home folder is "$HOMEFOLDER""
SEQS="data/seqs/"
ANALYSIS="analysis/"

cd ${HOMEFOLDER}${SEQS} # cd into the sequence folder

#### Read in sample list and make a bash array of the sample names (e.g. A1_S1)
find * -maxdepth 0 -name "*_L001_R1_001.fastq.gz" > samplelist.txt  # find all files ending with _L001_R1_001_fastq.gz
sample_info=samplelist.txt # put samplelist.txt into variable
sample_names=($(cut -f 1 "$sample_info" | uniq)) # convert variable to array this way
# echo ${sample_names[@]} # to echo all array elements
# echo ${sample_names[0]} # to echo first array element
echo "There are" ${#sample_names[@]} "samples that will be processed:  "${sample_names[@]} # echo number of elements in the array

# To run a loop interactively, select the entire loop and send to terminal.  Don't ctrl-Enter each line because this can send a command N times, as many lines as the command covers. So if the command wraps over 2 lines, ctrl-Entering each line sends the whole command twice, causing the command to be run twice per loop.

#### Use AdapterRemoval to trim Illumina sequencing adapters. (there are a few adapters still left in the raw data)
# brew install adapterremoval
# or download from https://github.com/MikkelSchubert/adapterremoval

# loop over all samples and remove any adapters
for sample in ${sample_names[@]}  # ${sample_names[@]} is the full bash array.  So loop over all samples
do
              sample_prefix="$( basename $sample "_L001_R1_001.fastq.gz")"
              echo ${sample_prefix}
              # AdapterRemoval --identify-adapters --file1 ${sample_prefix}_L001_R1_001.fastq.gz --file2 ${sample_prefix}_L001_R2_001.fastq.gz
              AdapterRemoval --file1 ${sample_prefix}_L001_R1_001.fastq.gz --file2 ${sample_prefix}_L001_R2_001.fastq.gz --basename Adaptermv${sample_prefix} --trimns
done

#### Sickle error correction strategies and identified quality trimming
# brew install sickle
# or download from https://github.com/najoshi/sickle
    # cd sickle-master
    # make
    # sudo cp sickle /usr/local/bin

# sickle -h

# loop over all samples and trim for quality
for sample in ${sample_names[@]}  # ${sample_names[@]} is the full bash array.  So loop over all samples
do
              sample_prefix="$( basename $sample "_L001_R1_001.fastq.gz")"  # e.g. A1_S1 is a sample_prefix
              echo ${sample_prefix} >> sickle.out
              echo ${sample_prefix}
              sickle pe -f Adaptermv${sample_prefix}.pair1.truncated -r Adaptermv${sample_prefix}.pair2.truncated -o sickle_${sample_prefix}_R1.fq -p sickle_${sample_prefix}_R2.fq -t sanger -s sickle_Single${sample_prefix}.fq >> sickle.out
done

##-s sickle_Single${sample_prefix}.fq:  keep single read （not paired but good sequences), will not be used in DAMe pipeline as we are using twin tags but can be used in genome research.

#### remove the AdapterRemoval files
rm Adaptermv*.*

#### SPAdes for error correction, via BayesHammer
#### to do:  test bfc
#### choosing not to use spades.py --meta because amplicon data
# download from http://cab.spbu.ru/software/spades/
    # sudo cp -r SPAdes-3.10.1 /usr/local/bin
# or
    # brew install homebrew/science/spades

# loop over all samples and apply BayesHammer error correction. Each file takes around 45 mins
for sample in ${sample_names[@]}  # ${sample_names[@]} is the full bash array.  So loop over all samples
do
              sample_prefix="$( basename $sample "_L001_R1_001.fastq.gz")"  # e.g. A1_S1 is a sample_prefix
              echo ${sample_prefix}
              date
              # spades.py --only-error-correction -1 sickle_${sample_prefix}_R1.fq -2 sickle_${sample_prefix}_R2.fq -s sickle_Single${sample_prefix}.fq -o SPAdes_hammer${sample_prefix}
              spades.py --only-error-correction -1 sickle_${sample_prefix}_R1.fq -2 sickle_${sample_prefix}_R2.fq -o SPAdes_hammer${sample_prefix} # don't denoise the unpaired reads
              date
done


#### remove the sickle files
rm sickle_*.fq

#### PandaSeq  use for read overlapping.  Each pair of fastq files takes ~ 15 secs
# download from https://github.com/neufeld/pandaseq/releases
    # PANDAseq-2.11.pkg
# or
    # brew install homebrew/science/pandaseq

for sample in ${sample_names[@]}  # ${sample_names[@]} is the full bash array.  So loop over all samples
do
    sample_prefix="$( basename $sample "_L001_R1_001.fastq.gz")"  # e.g. A1_S1 is a sample_prefix
    echo ${sample_prefix}
    date
    pandaseq -f SPAdes_hammer${sample_prefix}/corrected/sickle_${sample_prefix}_R1.00.0_0.cor.fastq.gz -r SPAdes_hammer${sample_prefix}/corrected/sickle_${sample_prefix}_R2.00.0_0.cor.fastq.gz -A simple_bayesian -d bfsrk -F -N -T 7 -g pandaseq_log_${sample_prefix}.txt -w sickle_cor_panda_${sample_prefix}.fastq
    date
                  # -A simple_bayesian # algorithm used in the original paper for pandaseq
                  # -F # output fastq file
                  # -N # eliminate all sequences with unknown nucleotides in the output
                  # -T 7 # use 7 threads
                  # -g pandaseq_log_${sample_prefix}.txt # output to this log file
                  # -d bfsrk # don't log this information (this set stops all logging i think)
done

# can replace the below with a parallel version. test command with --dryrun.  run for real without --dryrun
# However, this version uses the * wildcard, so it's not as precise as the loop
# find ./sickle_cor_panda_*.fastq -type f | parallel --dryrun gzip {}

# gzip the final fastq files
for sample in ${sample_names[@]}  # ${sample_names[@]} is the full bash array.  So loop over all samples
do
    sample_prefix="$( basename $sample "_L001_R1_001.fastq.gz")"  # e.g. A1_S1 is a sample_prefix
    echo ${sample_prefix}
    gzip sickle_cor_panda_${sample_prefix}.fastq
done

# remove SPAdes output
rm -rf SPAdes_hammer*/
# remove pandaseq_log_txt files
rm pandaseq_log_*.txt


#############################################################################################

#### DAMe - using PCR replicates and read numbers to filter out bad reads

# 1. place libraries in different folders
for sample in ${sample_names[@]}  # ${sample_names[@]} is the full bash array.  So loop over all samples
do
    sample_prefix="$( basename $sample "_L001_R1_001.fastq.gz")"  # e.g. A1_S1 is a sample_prefix
    echo ${sample_prefix}
    mkdir folder_${sample_prefix}
    mv sickle_cor_panda_${sample_prefix}.fastq.gz folder_${sample_prefix}/
done


# 2. SORT  Sort through each fastq file and determine how many of each tag pair is in each fastq file. Each fastq file takes < 1 min
for sample in ${sample_names[@]}  # ${sample_names[@]} is the full bash array.  So loop over all samples
do
              sample_prefix="$( basename $sample "_L001_R1_001.fastq.gz")"  # e.g. A1_S1 is a sample_prefix
              cd ${HOMEFOLDER}data/seqs # starting point
              echo ${sample_prefix}
              cd folder_${sample_prefix}
              python /usr/local/bin/DAMe/bin/sort.py -fq sickle_cor_panda_${sample_prefix}.fastq.gz -p ${HOMEFOLDER}data/Primers_COILeray.txt -t ${HOMEFOLDER}data/Tags_300test_COIA.txt
              # ls -lhS > sizesort_folder_${sample_prefix}.txt # quick check if the twin tag files are the largest. # sort by size.  The largest files should all be twin tag files (e.g. Tag1_Tag1.txt). Superseded by splitSummaryByPSInfo.py
done


# 3. Place PCR replicates of the same sample in the same folder and rename to pool{1,2,3}
# The three PCR replicate folders (on which sort.py was run) have to be in the same folder (e.g. 'folderA') and named 'pool1', 'pool2', and 'pool3'. No other pool folders

cd ${HOMEFOLDER}data/seqs

# original loop without if then for mkdir
# for sample in ${sample_names[@]}  # ${sample_names[@]} is the full bash array.  So loop over all samples
# do
#               sample_prefix="$( basename $sample "_L001_R1_001.fastq.gz")"
#               echo ${sample_prefix} | cut -c 1-2
#               # echo "${sample_prefix}" | cut -c n # selects out the nth character from sample_prefix
#               sample_prefix_prefix="$(echo "${sample_prefix}" | cut -c 1)"  # e.g. A is a sample_prefix_prefix
#               sample_prefix_pool="$(echo "${sample_prefix}" | cut -c 2)"  # e.g. 1 is a sample_prefix_pool
#               mkdir folder${sample_prefix_prefix}/  # make a folder with sample prefix
#                             # should replace above command with an if/then to mkdir folder only if it doesn't exist
#               mv folder_${sample_prefix}/ folder${sample_prefix_prefix}/pool${sample_prefix_pool}
# done


for sample in ${sample_names[@]}  # ${sample_names[@]} is the full bash array.  So loop over all samples
do
              sample_prefix="$( basename $sample "_L001_R1_001.fastq.gz")"
              echo ${sample_prefix} | cut -c 1-2
              # echo "${sample_prefix}" | cut -c n # selects out the nth character from sample_prefix
              sample_prefix_prefix="$(echo "${sample_prefix}" | cut -c 1)"  # e.g. A is a sample_prefix_prefix
              sample_prefix_pool="$(echo "${sample_prefix}" | cut -c 2)"  # e.g. 1 is a sample_prefix_pool
              if [ ! -d folder${sample_prefix_prefix} ] # if directory folder${sample_prefix_prefix} does not exist
              then
              	mkdir folder${sample_prefix_prefix}/  # make a folder with sample prefix
              fi
              mv folder_${sample_prefix}/ folder${sample_prefix_prefix}/pool${sample_prefix_pool}
done


# echo "${sample_prefix}" | cut -c 1-2 # selects out the first and second character from sample_prefix
# echo "${sample_prefix}" | grep -o '^\D' # selects first character, using grep

# 3.1 Sort SummaryCounts.txt to SummaryCountsSorted.txt (Bohmann code)

cd ${HOMEFOLDER}data/seqs

for sample in ${sample_names[@]}  # ${sample_names[@]} is the full bash array.  So loop over all samples
do
              sample_prefix="$( basename $sample "_L001_R1_001.fastq.gz")"
              cd ${HOMEFOLDER}data/seqs
              sample_prefix_prefix="$(echo "${sample_prefix}" | cut -c 1)"  # e.g. A is a sample_prefix_prefix
              sample_prefix_pool="$(echo "${sample_prefix}" | cut -c 2)"  # e.g. 1 is a sample_prefix_pool
              cd folder${sample_prefix_prefix}/pool${sample_prefix_pool}  # cd into each pool (e.g. folderA/pool1)
              head -1 SummaryCounts.txt > SummaryCounts_sorted_Folder${sample_prefix_prefix}_Pool${sample_prefix_pool}.txt
              tail -n +2 SummaryCounts.txt | sed "s/Tag//g" | sort -k1,1n -k2,2n | awk 'BEGIN{OFS="\t";}{$1="Tag"$1;$2="Tag"$2; print $0;}' >> SummaryCounts_sorted_Folder${sample_prefix_prefix}_Pool${sample_prefix_pool}.txt
              date
done

# 3.1.1 Move SummaryCounts_sorted*.txt files from inside pool folders into sample Folders
# find * -maxdepth 0 -name "*_L001_R1_001.fastq.gz" > samplelist.txt  # find all files ending with _L001_R1_001_fastq.gz
sample_libs=($(cat samplelist.txt | cut -c 1 | uniq))  # cut out all but the first letter of each filename and keep unique values, samplelist.txt should already exist
echo "There are" ${#sample_libs[@]} "samples that will be processed:  ${sample_libs[@]}." # echo number and name of elements in the array

# can replace the below with parallel syntax
# parallel echo "Moved SummaryCounts_sorted_Folder{1}_Pool{2}.txt" ::: ${sample_libs[@]} ::: `seq 1 ${POOLS}`

cd ${HOMEFOLDER}data/seqs
for sample in ${sample_libs[@]}  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
              for pool in `seq 1 ${POOLS}`
              do
                            echo "Moved SummaryCounts_sorted_Folder${sample}_Pool${pool}.txt"
                            mv folder${sample}/pool${pool}/SummaryCounts_sorted_Folder${sample}_Pool${pool}.txt folder${sample}
              done
done

# 3.2 Make tag combinations overview:  splitSummaryByPSInfo.py (Bohmann code)
# In the output file, the term "was used" means that the tag pair is in the PSInfo file (i.e. pair used in PCR)

cd ${HOMEFOLDER}data/seqs

for sample in ${sample_names[@]}  # ${sample_names[@]} is the full bash array.  So loop over all samples
do
              sample_prefix="$( basename $sample "_L001_R1_001.fastq.gz")"
              cd ${HOMEFOLDER}data/seqs
              sample_prefix_prefix="$(echo "${sample_prefix}" | cut -c 1)"  # e.g. A is a sample_prefix_prefix
              sample_prefix_pool="$(echo "${sample_prefix}" | cut -c 2)"  # e.g. 1 is a sample_prefix_pool
              # cd folder${sample_prefix_prefix}/pool${sample_prefix_pool}  # cd into each pool (e.g. folderA/pool1)
              cd folder${sample_prefix_prefix} # cd into each folder (e.g. folderA)
              python /usr/local/bin/DAMe/bin/splitSummaryByPSInfo.py -p ${HOMEFOLDER}data/PSinfo_300test_COI${sample_prefix_prefix}.txt -l ${sample_prefix_pool} -s pool${sample_prefix_pool}/SummaryCounts.txt -o splitSummaryByPSInfo_Folder${sample_prefix_prefix}_Pool${sample_prefix_pool}.txt
              date
done


# 3.3 Run scripts/heatmap.R on the different pools.

Rscript --vanilla --verbose ${HOMEFOLDER}scripts/heatmap.R

# then move heatmaps from inside pool folders into sample Folders
# should replace this with:  find | xargs mv
cd ${HOMEFOLDER}data/seqs
for sample in ${sample_libs[@]}  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
              for pool in `seq 1 ${POOLS}`
              do
                            mv folder${sample}/pool${pool}/heatmap_Folder${sample}_Pool${pool}.pdf folder${sample}
              done
done

# 3.9 Chimera checking [optional]  # I can't get this to run.  I get an error.  So i won't run it here.


# 4  Filter

# 4.1 RSI test PCR replicate similarity.  First filter.py at -y 1 -t 1, to keep all sequences, then run RSI.py on the pools
              ### This step is slow (~ 50 mins per library on my MacBookPro).  I run multiple jobs at once using parallel.

cd ${HOMEFOLDER}data/seqs
# python /usr/local/bin/DAMe/bin/filter.py -h
MINPCR_1=1 # min number of PCRs that a sequence has to appear in
MINREADS_1=1 # min number of copies per sequence per PCR
# confirm MINPCR and MINREADS values
echo "For RSI analysis, all sequences are kept: appear in just ${MINPCR_1} PCR, with just ${MINREADS_1} read per PCR."

#### This bit of code to make sample_libs[] exists above too (SummaryCountsSorted.txt).  Running here again too, in case, i'm running this loop independently
#### Read in sample list and make a bash array of the sample libraries (A, B, C, D, E, F)
# find * -maxdepth 0 -name "*_L001_R1_001.fastq.gz" > samplelist.txt  # find all files ending with _L001_R1_001_fastq.gz
sample_libs=($(cat samplelist.txt | cut -c 1 | uniq))  # cut out all but the first letter of each filename and keep unique values, samplelist.txt should already exist
# echo ${sample_libs[1]} # to echo first array element
echo "There are" ${#sample_libs[@]} "samples that will be processed:  ${sample_libs[@]}." # echo number and name of elements in the array

# GNU parallel install
              # brew install parallel
# https://unix.stackexchange.com/questions/294542/executing-commands-consequtively-on-multiple-folders
              # First, create a wrapper script that changes to the directory given in the first (and only) command-line argument, performs whatever setup/variable-initialisation/etc it needs, and then runs your 10 scripts in sequence with whatever args they need.
              #
              # For example, if each script processes all .jpg, .png, and .gif files in the directory:
              #
              # #! /bin/bash
              # # example-wrapper.sh
              #
              # cd "$1"
              #
              # script1 *.{jpg,png,gif}
              # script2 *.{jpg,png,gif}
              # script3 *.{jpg,png,gif}
              # script4 *.{jpg,png,gif}
              # script5 *.{jpg,png,gif}
              # script6 *.{jpg,png,gif}
              # script7 *.{jpg,png,gif}
              # script8 *.{jpg,png,gif}
              # script9 *.{jpg,png,gif}
              # script10 *.{jpg,png,gif}
              # Next, use find to pipe a list of directories into parallel.
              #
              # find /path/to/parent/ -mindepth 1 -type -d -print0 |
              #   parallel -0 -n 1 ./example-wrapper.sh
              # (the -mindepth 1 option in find excludes the top level directory, i.e. the parent directory itself)
              #
              # By default, parallel will run one instance (a "job") of ./example-wrapper.sh for each CPU core you have. Each instance will get ONE (-n 1) directory name. As soon as a job has finished, another is started (if there are any remaining jobs to run).
# https://stackoverflow.com/questions/11488184/how-can-i-run-a-list-of-commands-in-parallel
# https://stackoverflow.com/questions/11488184/how-can-i-run-a-list-of-commands-in-parallel
# alternative to parallel is rush, which takes multi-line commands, but it's too new yet:  https://github.com/shenwei356/rush/blob/master/README.md

cd ${HOMEFOLDER}data/seqs
# mkdir in parallel
parallel mkdir folder{1}/Filter_min${MINPCR_1}PCRs_min${MINREADS_1}copies_{1} ::: ${sample_libs[@]}
rm filter1_1_commands.txt # ensure no filter1_1_commands.txt file is present
# create a list of commands with the correct arguments
for sample in ${sample_libs[@]}  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
              echo "cd folder${sample}; \
              python /usr/local/bin/DAMe/bin/filter.py -psInfo ${HOMEFOLDER}data/PSinfo_300test_COI${sample}.txt -x ${PCRRXNS} -y ${MINPCR_1} -p ${POOLS} -t ${MINREADS_1} -l ${MINLEN} -o Filter_min${MINPCR_1}PCRs_min${MINREADS_1}copies_${sample}; \
              python /usr/local/bin/DAMe/bin/RSI.py --explicit -o RSI_output_${sample}.txt Filter_min${MINPCR_1}PCRs_min${MINREADS_1}copies_${sample}/Comparisons_${PCRRXNS}PCRs.txt" \
              >> filter1_1_commands.txt
done
# run parallel --dryrun to see the commands that will be run, without actually running them.
# parallel --dryrun --jobs 3 -k :::: filter1_1_commands.txt  # parallel :::: filter1_1_commands.txt means that the commands come from
parallel --jobs 3 -k :::: filter1_1_commands.txt  # parallel :::: filter1_1_commands.txt means that the commands come from filter1_1_commands.txt
rm filter1_1_commands.txt # ensure no filter1_1_commands.txt file is present

# non-parallel method.  analyses one folder at a time.
# cd ${HOMEFOLDER}data/seqs
# for sample in ${sample_libs[@]}  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
# do
#               cd ${HOMEFOLDER}data/seqs
#               cd folder${sample} # cd into folderA,B,C,D,E,F
#               mkdir Filter_min${MINPCR_1}PCRs_min${MINREADS_1}copies_${sample}
#               date
#
#               # filter.py
#               python /usr/local/bin/DAMe/bin/filter.py -psInfo ${HOMEFOLDER}data/PSinfo_300test_COI${sample}.txt -x ${PCRRXNS} -y ${MINPCR_1} -p ${POOLS} -t ${MINREADS_1} -l ${MINLEN} -o Filter_min${MINPCR_1}PCRs_min${MINREADS_1}copies_${sample}
#
#               # RSI.py
#               python /usr/local/bin/DAMe/bin/RSI.py --explicit -o RSI_output_${sample}.txt Filter_min${MINPCR_1}PCRs_min${MINREADS_1}copies_${sample}/Comparisons_${PCRRXNS}PCRs.txt
# done


# 4.2 plotLengthFreqMetrics_perSample.py and reads per sample per pool
# python /usr/local/bin/DAMe/bin/plotLengthFreqMetrics_perSample.py -h
cd ${HOMEFOLDER}data/seqs

for sample in ${sample_libs[@]}  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
              cd ${HOMEFOLDER}data/seqs
              cd folder${sample}/Filter_min${MINPCR_1}PCRs_min${MINREADS_1}copies_${sample} # cd into folderA,B,C,D,E,F/filter_output
              echo "Analysing: folder"${sample}/Filter_min${MINPCR_1}PCRs_min${MINREADS_1}copies_${sample}/
              # plotLengthFreqMetrics_perSample.py
              python /usr/local/bin/DAMe/bin/plotLengthFreqMetrics_perSample.py -f FilteredReads.fna -p ${HOMEFOLDER}data/PSinfo_300test_COI${sample}.txt -n 3
done

####################################################################################################
## After consideration of the negative controls, summary counts, and the heatmaps, choose thresholds for filtering
# The min PCR and copy number can be informed by looking at negative controls.
# After observing the negative controls, set -y and -t higher than the observed values in neg controls
####################################################################################################

# 4.3 Filter - Filtering reads with minPCR and minReads/PCR thresholds, min 300 bp length.
              ### This step is slow (~ 45 mins per library).

#### If I want to re-run filter.py with different thresholds, I set new values of MINPCR & MINREADS and start from here.
MINPCR=2 # min number of PCRs that a sequence has to appear in
MINREADS=4 # min number of copies per sequence per PCR
# confirm MINPCR and MINREADS values
echo "Each (unique) sequence must appear in at least ${MINPCR} PCRs, with at least ${MINREADS} reads per PCR."

cd ${HOMEFOLDER}data/seqs
# python /usr/local/bin/DAMe/bin/filter.py -h

#### Read in sample list and make a bash array of the sample libraries (A, B, C, D, E, F)
# find * -maxdepth 0 -name "*_L001_R1_001.fastq.gz" > samplelist.txt  # find all files ending with _L001_R1_001_fastq.gz
sample_libs=($(cat samplelist.txt | cut -c 1 | uniq))  # cut out all but the first letter of each filename and keep unique values, samplelist.txt should already exist
# echo ${sample_libs[1]} # to echo first array element
echo "There are" ${#sample_libs[@]} "samples that will be processed:  ${sample_libs[@]}." # echo number and name of elements in the array

for sample in ${sample_libs[@]}  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
              cd ${HOMEFOLDER}data/seqs
              cd folder${sample} # cd into folderA,B,C,D,E,F
              mkdir Filter_min${MINPCR}PCRs_min${MINREADS}copies_${sample}
              date
              python /usr/local/bin/DAMe/bin/filter.py -psInfo ${HOMEFOLDER}data/PSinfo_300test_COI${sample}.txt -x ${PCRRXNS} -y ${MINPCR} -p ${POOLS} -t ${MINREADS} -l ${MINLEN} -o Filter_min${MINPCR}PCRs_min${MINREADS}copies_${sample}
done
              ## python /usr/local/bin/DAMe/bin/filter.py -h
               #  -x X                  Number of PCR rxns performed per sample
               #  -y Y                  Number of PCR rxns in which the sequence has to be
               #                        present
               #  -p P                  The number of pools in which the samples were divided
               #                        for sequencing (in case of tag combinations repetition
               #                        due to processing so many samples) [default 1] NOTE:
               #                        If using pools, each fastq must be in a folder called
               #                        pool#, in which the sort.py was run for each pool
               #                        inside the corresponding folder, and this program
               #                        chimeraCheck.py is run in the parent directory of the
               #                        pools directories
               #  -t T                  Number of times a unique sequence has to be present in a PCR rxn
               #  -l L                  Minimum sequence length
               # --chimeraChecked      Use this parameter if you have performed a chimera
               #                          check on the sorted collapsed sequence files [default
               #                          not set]
               #  -o OUT, --outDir OUT  Output directory

# 4.4 plotLengthFreqMetrics_perSample.py and reads per sample per pool
# python /usr/local/bin/DAMe/bin/plotLengthFreqMetrics_perSample.py -h

cd ${HOMEFOLDER}data/seqs

for sample in ${sample_libs[@]}  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
              cd ${HOMEFOLDER}data/seqs
              cd folder${sample}/Filter_min${MINPCR}PCRs_min${MINREADS}copies_${sample} # cd into folderA,B,C,D,E,F/filter_output
              echo "Analysing: folder"${sample}/Filter_min${MINPCR}PCRs_min${MINREADS}copies_${sample}/
              # plotLengthFreqMetrics_perSample.py
              python /usr/local/bin/DAMe/bin/plotLengthFreqMetrics_perSample.py -f FilteredReads.fna -p ${HOMEFOLDER}data/PSinfo_300test_COI${sample}.txt -n 3
done


# 5. Prepare the FilteredReads.fna file for Sumaclust clustering. changes the header lines on the fasta file
# install sumatra
              # cd ~/src
              # wget https://git.metabarcoding.org/obitools/sumatra/uploads/251020bbbd6c6595cb9fce6077e29952/sumatra_v1.0.20.tar.gz
              # tar -zxvf sumatra_v1.0.20.tar.gz
              # cd sumatra_v1.0.20/
              # make CC=clang # disables OpenMP, which isn't on macOS
              # mv sumatra /usr/local/bin

# 5.1 assessing clustering parameters
# This takes quite a long time, only feasible on FilteredReads.fna that have been filtered to a few MB
cd ${HOMEFOLDER}data/seqs/
# MINPCR=2 # these commands are to make it possible to process multiple filter.py outputs, which are saved in different Filter_min folders
# MINREADS=4
echo "Analysing filter.py output where each (unique) sequence appeared in ≥ ${MINPCR} PCRs, with ≥ ${MINREADS} reads per PCR."

for sample in ${sample_libs[@]}  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
              cd ${HOMEFOLDER}data/seqs/folder${sample}/Filter_min${MINPCR}PCRs_min${MINREADS}copies_${sample}
              python /usr/local/bin/DAMe/bin/assessClusteringParameters.py -i FilteredReads.fna -mint 0.8 -minR 0.6 -step 0.05 -t 4 -o 16sclusterassess_mint08_minR06_step005_Filter_min${MINPCR}PCRs_min${MINREADS}copies_${sample}.pdf
done
              # python /usr/local/bin/DAMe/bin/assessClusteringParameters.py -h
              # Generate plots to explore the parameter space for OTU clustering parameters.
              # optional arguments:
              #   -h, --help            show this help message and exit
              #   -i InputFasta, --inFasta InputFasta
              #                         Input fasta file for cluster.
              #   -mint minT, --minIdentityCutoff minT
              #                         Minimum identity cutoff for sumatra/sumaclust
              #   -minR minR, --minAbundanceRatio minR
              #                         Minimum abundance ratio for sumatra/sumaclust
              #   -step stepSize, --stepSize stepSize
              #                         Step size for t and R
              #   -t threads, --threads threads
              #                         Number of threads to use
              #   -o outfile, --out outfile
              #                         Output file pdf


# 5.2 python /usr/local/bin/DAMe/bin/convertToUSearch.py -h

# MINPCR=2 # these commands are to make it possible to process multiple filter.py outputs, which are saved in different Filter_min folders
# MINREADS=3

for sample in ${sample_libs[@]}  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
              cd ${HOMEFOLDER}data/seqs/folder${sample}/Filter_min${MINPCR}PCRs_min${MINREADS}copies_${sample}
              python /usr/local/bin/DAMe/bin/convertToUSearch.py -i FilteredReads.fna -lmin 300 -lmax 320
              python /usr/local/bin/DAMe/bin/convertToUSearch.py -i FilteredReads.fna -lmin 300 -lmax 320 -u
done

#### If I instead run with --usearch, I can then use a usearch pipeline to search for chimeras, cluster OTUs (ZOTUs in usearch parlance), and generate OTU tables.  I can't get this to work
              # python /usr/local/bin/DAMe/bin/convertToUSearch.py -i FilteredReads.fna -lmin 300 -lmax 320 --usearch
              # cd ${HOMEFOLDER}data/seqs/folderA
              # date; usearch10 -uchime3_denovo FilteredReads.forusearch.fas -uchimeout out.txt -chimeras ch.fa -nonchimeras FilteredReads.forusearch.nochm.fna; date


# 6. Sumaclust clustering and convert Sumaclust output to table format

# sumaclust download from: https://git.metabarcoding.org/obitools/sumaclust/wikis/home
# wget https://git.metabarcoding.org/obitools/sumaclust/uploads/69f757c42f2cd45212c587e87c75a00f/sumaclust_v1.0.20.tar.gz
# tar -zxvf sumaclust_v1.0.20.tar.gz
# cd sumaclust_v1.0.20/
# make CC=clang # disables OpenMP, which isn't on macOS
# mv sumaclust /usr/local/bin/


# ~/src/sumaclust_v1.0.20/sumaclust -h
# python /usr/local/bin/DAMe/bin/tabulateSumaclust.py -h

# 96% sumaclust
echo ${SUMASIM} # confirm that there is a similarity value chosen
SUMASIM=96 # if there is no SUMASIM value
echo ${SUMASIM} # confirm that there is a similarity value chosen
cd ${HOMEFOLDER}
for sample in ${sample_libs[@]}  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
              cd ${HOMEFOLDER}data/seqs/folder${sample}/Filter_min${MINPCR}PCRs_min${MINREADS}copies_${sample}
              sumaclust -t .${SUMASIM} -e FilteredReads.forsumaclust.fna > OTUs_${SUMASIM}_sumaclust.fna
              python /usr/local/bin/DAMe/bin/tabulateSumaclust.py -i OTUs_${SUMASIM}_sumaclust.fna -o table_300test_${sample}_${SUMASIM}.txt -blast
done

# 97% sumaclust
# MINPCR=2 # these commands are to make it possible to prepare multiple filter.py outputs
# MINREADS=3
echo ${SUMASIM} # confirm the similarity value
SUMASIM=97 # if there is no SUMASIM value
echo ${SUMASIM} # confirm the similarity value
cd ${HOMEFOLDER}
for sample in ${sample_libs[@]}  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
              cd ${HOMEFOLDER}data/seqs/folder${sample}/Filter_min${MINPCR}PCRs_min${MINREADS}copies_${sample}
              sumaclust -t .${SUMASIM} -e FilteredReads.forsumaclust.fna > OTUs_${SUMASIM}_sumaclust.fna
              python /usr/local/bin/DAMe/bin/tabulateSumaclust.py -i OTUs_${SUMASIM}_sumaclust.fna -o table_300test_${sample}_${SUMASIM}.txt -blast
done



# 7. Move sumaclust results to ${HOMEFOLDER}/analysis

# MINPCR=2; MINREADS=3 # if not running as bash

for sample in ${sample_libs[@]}  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
              cd ${HOMEFOLDER}${ANALYSIS}
              # make folder to hold results
              if [ ! -d OTU_transient_results ] # if directory OTU_transient_results does not exist
              then
              	mkdir OTU_transient_results/
                            mkdir OTU_transient_results/OTU_tables/
              fi
              mv ${HOMEFOLDER}${SEQS}folder${sample}/Filter_min${MINPCR}PCRs_min${MINREADS}copies_${sample} ${HOMEFOLDER}${ANALYSIS}OTU_transient_results/Filter_min${MINPCR}PCRs_min${MINREADS}copies_${sample}
              cd ${HOMEFOLDER}${ANALYSIS}OTU_transient_results/Filter_min${MINPCR}PCRs_min${MINREADS}copies_${sample}
              rm Comparisons_${PCRRXNS}PCRs.fasta # large file that doesn't need to be kept
              rm Comparisons_${PCRRXNS}PCRs.txt # large file that doesn't need to be kept
              cp table_*_${sample}_*.txt ${HOMEFOLDER}${ANALYSIS}OTU_transient_results/OTU_tables/ # copy OTU tables to OTU_tables folder
done

# change name of OTU_transient_results folder to include timestamp
mv ${HOMEFOLDER}${ANALYSIS}OTU_transient_results/ ${HOMEFOLDER}${ANALYSIS}OTU_tables+seqs_$(date +%F_time-%H%M)/

# This is where i should run usearch9 -uchime2_denovo
#### End script here.  Everything below is not to be run
exit




#### ***** Using CROP makes 346 97% OTUs, versus 391 sumaclust 97% OTUs
# cd folderA/Filter_min2PCRs_min2copies_A/
# ~/src/CROP/CROP -i FilteredReads.forsumaclust.fna -o OTUs_CROP_97.out -s -e 2000 -b 185 -z 300 -r 0
              # Parms are:
              # -r = "This parameter controls the size of the clusters to be considered as “rare”.  –r 0 will yield the best accuracy
              # -z < 150,000/avg_seq_length = 150,000/320 bp = 469. let zl < 150000, so set z=300;
              # -b = block size = number_of_sequences/50 = 9227/50 = 185;  300Test/folderA/Filter_min2PCRs_min2copies_A/FilteredReads.forsumaclust.fna
              # -e = 10*block size = 1850.  2000 (default)
              # -s # 97% similarity

#  HOWEVER, CROP takes 5 minutes, and sumaclust takes 10 secs. Also, using CROP needs a need script to convert the CROP OTU map (OTUs_CROP_97.out.cluster.list) and the OTU representative seq list (OTUs_CROP_97.out.cluster.fasta) into an OTU table that preserves the number of original Illumina reads per OTU-sample combination (e.g. table_300test_A_97.txt, which is post sumaclust + tabulateSumaclust.py + OTUs_97_sumaclust.fna)

##### MTB and PC Sanger sequences:  MTBallSpeciesSeqs_20170719.fasta

cd ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/
~/src/sumaclust_v1.0.20/sumaclust -h

~/src/sumaclust_v1.0.20/sumaclust -t 0.98 -e MTBallSpeciesSeqs_20170719.fasta > MTB_OTUs_98_sumaclust.fas
# 261 clusters
~/src/sumaclust_v1.0.20/sumaclust -t 0.97 -e MTBallSpeciesSeqs_20170719.fasta > MTB_OTUs_97_sumaclust.fas
# 254 clusters
~/src/sumaclust_v1.0.20/sumaclust -t 0.96 -e MTBallSpeciesSeqs_20170719.fasta > MTB_OTUs_96_sumaclust.fas
# 252 clusters
~/src/CROP/CROP -i MTBallSpeciesSeqs_20170719.fasta -o MTB_OTUs_CROP_97.out -s -e 60 -b 6 -z 230 -r 0
# 231 clusters

# Parms are:
# -r = "This parameter controls the size of the clusters to be considered as “rare”.  –r 0 will yield the best accuracy
# -z < 150,000/avg_seq_length = 150,000/623 bp = 241. Let zl < 150,000, so set z=230;
# -b = block size = number_of_sequences/50 = 280/50 = 5.6, so set to b=6
# -e = 10*block size = 60.  2000 (default)
# -s # 97% similarity


#########################
7. Choose final filter.py thresholds

cd ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/
mkdir otutables/
cp folder{A,B,C,D,E,F}/Filter_min2PCRs_min2copies_{A,B,C,D,E,F}/table_300test_{A,B,C,D,E,F}_96.txt ./otutables
cp folder{A,B,C,D,E,F}/Filter_min2PCRs_min2copies_{A,B,C,D,E,F}/table_300test_{A,B,C,D,E,F}_97.txt ./otutables
# insert into an excel table and check for mix-ups, extraction blanks,  positive control, and samples. decide on new thresholds for filtering
# consider using R to read the tables and write to an excel file

##filtering with new thresholds
mkdir Filter_min2PCRs_min10copies_A
date; python /usr/local/bin/DAMe/bin/filter.py -psInfo PSinfo_300test_COIA.txt -x 3 -y 2 -p 3 -t 10 -l 300 -o Filter_min2PCRs_min20copies_A ;date
#2017年 4月14日 星期五 22时02分06秒 CST
#2017年 4月14日 星期五 22时46分19秒 CST
python /usr/local/bin/DAMe/bin/convertToUSearch.py -h
cd Filter_min2PCRs_min10copies_A
python /usr/local/bin/DAMe/bin/convertToUSearch.py -i FilteredReads.fna -lmin 300 -lmax 320
sumaclust -t 0.96 -e FilteredReads.forsumaclust.fna > OTUs_96_sumaclust.fna
python /usr/local/bin/DAMe/bin/tabulateSumaclust.py -i OTUs_96_sumaclust.fna -o table_300test_A_96.txt -blast


mkdir Filter_min2PCRs_min3copies_A
date; python /usr/local/bin/DAMe/bin/filter.py -psInfo PSinfo_300test_COIA.txt -x 3 -y 2 -p 3 -t 3 -l 300 -o Filter_min2PCRs_min3copies_A ;date
#2017年 4月14日 星期五 22时11分11秒 CST
#2017年 4月14日 星期五 22时54分34秒 CST
cd Filter_min2PCRs_min3copies_A
python /usr/local/bin/DAMe/bin/convertToUSearch.py -i FilteredReads.fna -lmin 300 -lmax 320
sumaclust -t 0.96 -e FilteredReads.forsumaclust.fna > OTUs_96_sumaclust.fna
python /usr/local/bin/DAMe/bin/tabulateSumaclust.py -i OTUs_96_sumaclust.fna -o table_300test_A_96.txt -blast


mkdir Filter_min3PCRs_min3copies_A
date; python /usr/local/bin/DAMe/bin/filter.py -psInfo PSinfo_300test_COIA.txt -x 3 -y 3 -p 3 -t 3 -l 300 -o Filter_min3PCRs_min3copies_A ;date
#2017年 4月14日 星期五 23时05分11秒 CST
#2017年 4月14日 星期五 23时46分21秒 CST
cd Filter_min3PCRs_min3copies_A
python /usr/local/bin/DAMe/bin/convertToUSearch.py -i FilteredReads.fna -lmin 300 -lmax 320
sumaclust -t 0.96 -e FilteredReads.forsumaclust.fna > OTUs_96_sumaclust.fna  ###276  OTU
python /usr/local/bin/DAMe/bin/tabulateSumaclust.py -i OTUs_96_sumaclust.fna -o table_300test_A_96.txt -blast


mkdir Filter_min3PCRs_min2copies_A
date; python /usr/local/bin/DAMe/bin/filter.py -psInfo PSinfo_300test_COIA.txt -x 3 -y 3 -p 3 -t 2 -l 300 -o Filter_min3PCRs_min2copies_A ;date
#2017年 5月 5日 星期五 11时10分08秒 CST
#2017年 5月 5日 星期五 11时55分52秒 CST
cd Filter_min3PCRs_min2copies_A
python /usr/local/bin/DAMe/bin/convertToUSearch.py -i FilteredReads.fna -lmin 300 -lmax 320
sumaclust -t 0.96 -e FilteredReads.forsumaclust.fna > OTUs_96_sumaclust.fna  ###301 OTU
python /usr/local/bin/DAMe/bin/tabulateSumaclust.py -i OTUs_96_sumaclust.fna -o table_300test_A_96.txt -blast

mkdir Filter_min3PCRs_min5copies_A
date; python /usr/local/bin/DAMe/bin/filter.py -psInfo PSinfo_300test_COIA.txt -x 3 -y 3 -p 3 -t 5 -l 300 -o Filter_min3PCRs_min5copies_A ;date
#2017年 4月29日 星期六 08时54分05秒 CST
#2017年 4月29日 星期六 09时56分39秒 CST
cd Filter_min3PCRs_min5copies_A
python /usr/local/bin/DAMe/bin/convertToUSearch.py -i FilteredReads.fna -lmin 300 -lmax 320
sumaclust -t 0.96 -e FilteredReads.forsumaclust.fna > OTUs_96_sumaclust.fna  ###260 OTU
python /usr/local/bin/DAMe/bin/tabulateSumaclust.py -i OTUs_96_sumaclust.fna -o table_300test_A_96.txt -blast


cd B
mkdir Filter_min3PCRs_min2copies_B
date; python /usr/local/bin/DAMe/bin/filter.py -psInfo PSinfo_300test_COIB.txt -x 3 -y 3 -p 3 -t 2 -l 300 -o Filter_min3PCRs_min2copies_B ;date
#2017年 4月19日 星期三 11时01分37秒 CST
#2017年 4月19日 星期三 11时39分20秒 CST
cd Filter_min3PCRs_min2copies_B
python /usr/local/bin/DAMe/bin/convertToUSearch.py -i FilteredReads.fna -lmin 300 -lmax 320
sumaclust -t 0.96 -e FilteredReads.forsumaclust.fna > OTUs_96_sumaclust.fna  ###305 OTU
python /usr/local/bin/DAMe/bin/tabulateSumaclust.py -i OTUs_96_sumaclust.fna -o table_300test_A_96.txt -blast

mkdir Filter_min3PCRs_min1copies_B
date; python /usr/local/bin/DAMe/bin/filter.py -psInfo PSinfo_300test_COIB.txt -x 3 -y 3 -p 3 -t 2 -l 300 -o Filter_min3PCRs_min1copies_B ;date
#2017年 4月19日 星期三 13时03分59秒 CST
#2017年 4月19日 星期三 13时44分00秒 CST
cd Filter_min3PCRs_min1copies_B
python /usr/local/bin/DAMe/bin/convertToUSearch.py -i FilteredReads.fna -lmin 300 -lmax 320
sumaclust -t 0.96 -e FilteredReads.forsumaclust.fna > OTUs_96_sumaclust.fna  ###305 OTU
python /usr/local/bin/DAMe/bin/tabulateSumaclust.py -i OTUs_96_sumaclust.fna -o table_300test_A_96.txt -blast

mkdir Filter_min2PCRs_min2copies_B
date; python /usr/local/bin/DAMe/bin/filter.py -psInfo PSinfo_300test_COIB.txt -x 3 -y 3 -p 3 -t 2 -l 300 -o Filter_min2PCRs_min2copies_B ;date
#2017年 4月19日 星期三 13时57分54秒 CST
#2017年 4月19日 星期三 14时38分28秒 CST
cd Filter_min2PCRs_min2copies_B
python /usr/local/bin/DAMe/bin/convertToUSearch.py -i FilteredReads.fna -lmin 300 -lmax 320
sumaclust -t 0.96 -e FilteredReads.forsumaclust.fna > OTUs_96_sumaclust.fna  ###305  OTU
python /usr/local/bin/DAMe/bin/tabulateSumaclust.py -i OTUs_96_sumaclust.fna -o table_300test_A_96.txt -blast

mkdir Filter_min3PCRs_min5copies_B
date; python /usr/local/bin/DAMe/bin/filter.py -psInfo PSinfo_300test_COIB.txt -x 3 -y 3 -p 3 -t 5 -l 300 -o Filter_min3PCRs_min5copies_B ;date
#2017年 4月29日 星期六 08时55分52秒 CST
#2017年 4月29日 星期六 09时53分55秒 CST
cd Filter_min3PCRs_min5copies_B
python /usr/local/bin/DAMe/bin/convertToUSearch.py -i FilteredReads.fna -lmin 300 -lmax 320
sumaclust -t 0.96 -e FilteredReads.forsumaclust.fna > OTUs_96_sumaclust.fna  ###266  OTU
python /usr/local/bin/DAMe/bin/tabulateSumaclust.py -i OTUs_96_sumaclust.fna -o table_300test_A_96.txt -blast


cd C
mkdir Filter_min3PCRs_min2copies_C
date; python /usr/local/bin/DAMe/bin/filter.py -psInfo PSinfo_300test_COIC.txt -x 3 -y 3 -p 3 -t 2 -l 300 -o Filter_min3PCRs_min2copies_C ;date
#2017年 4月20日 星期四 11时55分56秒 CST
#2017年 4月20日 星期四 12时39分56秒 CST
cd Filter_min2PCRs_min2copies_C
python /usr/local/bin/DAMe/bin/convertToUSearch.py -i FilteredReads.fna -lmin 300 -lmax 320
sumaclust -t 0.96 -e FilteredReads.forsumaclust.fna > OTUs_96_sumaclust.fna  ###281  OTU
python /usr/local/bin/DAMe/bin/tabulateSumaclust.py -i OTUs_96_sumaclust.fna -o table_300test_A_96.txt -blast

mkdir Filter_min3PCRs_min5copies_C
date; python /usr/local/bin/DAMe/bin/filter.py -psInfo PSinfo_300test_COIC.txt -x 3 -y 3 -p 3 -t 5 -l 300 -o Filter_min3PCRs_min5copies_C ;date
#2017年 4月29日 星期六 09时37分23秒 CST
#2017年 4月29日 星期六 10时30分18秒 CST
cd Filter_min3PCRs_min5copies_C
python /usr/local/bin/DAMe/bin/convertToUSearch.py -i FilteredReads.fna -lmin 300 -lmax 320
sumaclust -t 0.96 -e FilteredReads.forsumaclust.fna > OTUs_96_sumaclust.fna  ###248  OTU
python /usr/local/bin/DAMe/bin/tabulateSumaclust.py -i OTUs_96_sumaclust.fna -o table_300test_A_96.txt -blast


cd D
mkdir Filter_min3PCRs_min2copies_D
date; python /usr/local/bin/DAMe/bin/filter.py -psInfo PSinfo_300test_COID.txt -x 3 -y 3 -p 3 -t 2 -l 300 -o Filter_min3PCRs_min2copies_D ;date
2017年 4月21日 星期五 17时59分16秒 CST
2017年 4月21日 星期五 18时49分42秒 CST
cd Filter_min3PCRs_min2copies_D
python /usr/local/bin/DAMe/bin/convertToUSearch.py -i FilteredReads.fna -lmin 300 -lmax 320
sumaclust -t 0.96 -e FilteredReads.forsumaclust.fna > OTUs_96_sumaclust.fna  ###273  OTU
python /usr/local/bin/DAMe/bin/tabulateSumaclust.py -i OTUs_96_sumaclust.fna -o table_300test_A_96.txt -blast

mkdir Filter_min3PCRs_min5copies_D
date; python /usr/local/bin/DAMe/bin/filter.py -psInfo PSinfo_300test_COID.txt -x 3 -y 3 -p 3 -t 5 -l 300 -o Filter_min3PCRs_min5copies_D ;date
#2017年 4月29日 星期六 10时19分06秒 CST
#2017年 4月29日 星期六 11时12分21秒 CST
cd Filter_min3PCRs_min5copies_D
python /usr/local/bin/DAMe/bin/convertToUSearch.py -i FilteredReads.fna -lmin 300 -lmax 320

cd E
mkdir Filter_min3PCRs_min2copies_E
date; python /usr/local/bin/DAMe/bin/filter.py -psInfo PSinfo_300test_COIE.txt -x 3 -y 3 -p 3 -t 2 -l 300 -o Filter_min3PCRs_min2copies_E ;date
#2017年 4月21日 星期五 20时16分16秒 CST
#2017年 4月21日 星期五 20时51分48秒 CST
cd Filter_min3PCRs_min2copies_E
python /usr/local/bin/DAMe/bin/convertToUSearch.py -i FilteredReads.fna -lmin 300 -lmax 320
sumaclust -t 0.96 -e FilteredReads.forsumaclust.fna > OTUs_96_sumaclust.fna  ###302  OTU
python /usr/local/bin/DAMe/bin/tabulateSumaclust.py -i OTUs_96_sumaclust.fna -o table_300test_A_96.txt -blast

mkdir Filter_min3PCRs_min5copies_E
date; python /usr/local/bin/DAMe/bin/filter.py -psInfo PSinfo_300test_COIE.txt -x 3 -y 3 -p 3 -t 5 -l 300 -o Filter_min3PCRs_min5copies_E ;date


cd F
mkdir Filter_min3PCRs_min2copies_F
date; python /usr/local/bin/DAMe/bin/filter.py -psInfo PSinfo_300test_COIF.txt -x 3 -y 3 -p 3 -t 2 -l 300 -o Filter_min3PCRs_min2copies_F ;date
#2017年 4月22日 星期六 13时51分56秒 CST
#2017年 4月22日 星期六 14时15分57秒 CST
cd Filter_min3PCRs_min2copies_F
python /usr/local/bin/DAMe/bin/convertToUSearch.py -i FilteredReads.fna -lmin 300 -lmax 320
sumaclust -t 0.96 -e FilteredReads.forsumaclust.fna > OTUs_96_sumaclust.fna  ###289  OTU
python /usr/local/bin/DAMe/bin/tabulateSumaclust.py -i OTUs_96_sumaclust.fna -o table_300test_A_96.txt -blast

mkdir Filter_min3PCRs_min5copies_F
date; python /usr/local/bin/DAMe/bin/filter.py -psInfo PSinfo_300test_COIF.txt -x 3 -y 3 -p 3 -t 5 -l 300 -o Filter_min3PCRs_min5copies_F ;date
#2017年 4月29日 星期六 08时52分35秒 CST
#2017年 4月29日 星期六 09时27分47秒 CST
cd Filter_min3PCRs_min5copies_F
python /usr/local/bin/DAMe/bin/convertToUSearch.py -i FilteredReads.fna -lmin 300 -lmax 320
sumaclust -t 0.96 -e FilteredReads.forsumaclust.fna > OTUs_96_sumaclust.fna  ###248  OTU
python /usr/local/bin/DAMe/bin/tabulateSumaclust.py -i OTUs_96_sumaclust.fna -o table_300test_A_96.txt -blast





#####建立本地blast 数据库？？蔡望做
makeblastdb -in MTBallSpeciesSeqs_noYWK.fasta -dbtype nucl

usearch9 -sortbysize unidentified_OTU_B.fasta -fastaout unidentified_OTU_B_sorted.fasta

usearch9 -uchime2_ref unidentified_OTU_B.fasta -db MTBallSpeciesSeqs_noYWK.fasta -uchimeout out.txt -strand plus -mode sensitive

usearch9 -uchime2_denovo unidentified_OTU_B_sorted.fasta -chimeras denovo_chimeras9.fasta



############################################比较不同的Tag差异########################################
###Use the output “sickle_cor_panda_XX.fna”after "AdapterRemoval+Sickle+SPades+Pandaseq"###
###on QIIME
source /macqiime/configs/bash_profile.txt		### for using Macqiime

###PCR_A1
validate_mapping_file.py -m Mappingfile_A1.txt -o check_map/

#convert fastq file to fasta file.
convert_fastaqual_fastq.py -f sickle_cor_panda_A1.fastq -c fastq_to_fastaqual

date; split_libraries.py -m Mappingfile_A1.txt -f sickle_cor_panda_A1.fna -q sickle_cor_panda_A1.qual -o 2split_A1/ -l 300 -L 450 -b hamming_8 -r -t -d -z truncate_remove -s 20 -M 2 --reverse_primer_mismatches 2;date


###use seqtk looking for reverse matches
seqtk seq -r sickle_cor_panda_A1.fastq >sickle_cor_panda_A1_rc.fastq
convert_fastaqual_fastq.py -f sickle_cor_panda_A1_rc.fastq -c fastq_to_fastaqual
date; split_libraries.py -m Mappingfile_A1.txt -f sickle_cor_panda_A1_rc.fna -q sickle_cor_panda_A1_rc.qual -o 2split_A1_rc/ -l 300 -L 450 -b hamming_8 -r -t -d -z truncate_remove -s 20 -M 2 --reverse_primer_mismatches 2 -n 2464;date


#########################################比较不同的Tag差异，用Wang Xiaoyang的perl程序####################
#程序名称：300Test_extractAllreadsfromSort.pl
#Author: Wang Xiaoyang
#Aim: put all seqs from different pools into one single fasta file.

perl 300Test_extractAllreadsfromSort.pl -id AF_Folder/ -o AF_allseqs.fasta

#########################################uchime to remove chimera################################

date; usearch9 -fastx_uniques AF_allseqs.fasta -fastaout AF_allseqs_uniques.fasta -sizeout -relabel unique -uc unique.uc; date
###1390945 uniques.  1205789 singleton reads (86.7%)

###     just used unique.uc file; put uc.files into one folder--uc
perl Pick_cluster_num_V3.pl -id uc/ -o usearch9_map.txt		## added function for getting size of each unique seq

### pick rep seqs
date; pick_rep_set.py -i usearch9_map.txt -f AF_allseqs.fasta -o usearch_dereped.fasta -m most_abundant; date
2017年 7月12日 星期三 11时56分16秒 CST
2017年 7月12日 星期三 12时03分53秒 CST

### deChimera

#uchime2 with unique sequences
less usearch_dereped.fasta | tr ' ' ';' > usearch_dereped9.fasta   ###replace' 'with';'
usearch9 -sortbysize usearch_dereped9.fasta -fastaout usearch_dereped9_sorted.fasta
# if run on mac
date; usearch9 -uchime2_denovo usearch_dereped9_sorted.fasta -chimeras denovo_chimeras9.fasta -nonchimeras denovo_nonchimeras9.fasta; date  # if run on mac.

# This does not work yet.
# date; usearch10 -uchime3_denovo denoised.fa -uchimeout out.txt -chimeras ch.fa -nonchimeras nonch.fa; date


# if run on server
ssh wangxy@10.0.16.80
kiz650223
cd Yahan
screen
date; ~/scripts/usearch/usearch9.2.64_i86linux32 -uchime2_denovo usearch_dereped9_sorted.fasta -chimeras denovo_chimeras9.fasta  -nonchimeras denovo_nonchimeras9.fasta; date
# ctrl-a d  # type this into the terminal
screen -r 15693

# uchime2 without unique sequences
# create a fasta file without singleton reads (unique reads)
~/scripts/usearch/usearch9.2.64_i86linux32 -sortbysize usearch_dereped9_sorted.fasta -fastaout usearch_dereped9_sorted_nouniques.fasta -minsize 2
# usearch -fastx_uniques input.fasta -fastaout uniques.fasta -minuniquesize 2 -sizeout # alternative method
# ./usearch10 -fastx_info usearch_dereped9_sorted_nouniques.fasta
# File size 77.6M, 215.6k seqs, 67.4M letters
# Lengths min 15, low 312, med 313, hi 313, max 442
# Letter freqs T 38.7%, A 28.9%, C 18.4%, G 14.0%
# 0% masked (lower-case)

screen
date; ~/scripts/usearch/usearch9.2.64_i86linux32 -uchime2_denovo usearch_dereped9_sorted_nouniques.fasta -chimeras denovo_nouniques_chimeras9.fasta  -nonchimeras denovo_nouniques_nonchimeras9.fasta; date
# ctrl-a d  # type this into the terminal
screen -r 15791
