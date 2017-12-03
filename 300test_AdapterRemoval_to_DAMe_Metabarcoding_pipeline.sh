#!/bin/bash
set -e
set -u
set -o pipefail
##########################################################################################################
##########################################################################################################
# shell script for metabarcoding:  Schirmer et al (2015, NAR) + Zepeda-Mendoza et al (2016) DAMe pipeline
##########################################################################################################
###########################################################################################################

# To do
              # try bfc instead of spades.py --error_correction-only
              # understand output of assessClusteringParameters.py


# Usage: bash 300test_AdaptorRMV+Panda+DAMe.sh SUMASIM
# e.g. 300test_AdaptorRMV+Panda+DAMe.sh 2 3 97 # OTU must appear in at least 2 PCRs and 3 reads per PCR, and will be clustered at 97%

# PIPESTART=$(date)

# run script from ~/PATH_TO/300Test/scripts
# fastq files in 300Test/data/seqs/
# DAMe files in 300Test/data/
# analysis outputs in 300Test/analysis/
# this script in 300Test/scripts/

# set variables
# SUMASIM=$1 # sumaclust similarity in percentage (0-100, e.g. 97)
SUMASIM=97 # sumaclust similarity in percentage (0-100, e.g. 97)
MINPCR=2 # min contig length after assembly if i'm not running as a bash script
MINREADS=4 # min contig length after assembly if i'm not running as a bash script
MINPCR_1=1 # for RSI analysis and to run without filtering by min PCR number or min copy number per pcr
MINREADS_1=1 # for RSI analysis and to run without filtering by min PCR number or min copy number per pcr
MINLEN=300
MAXLEN=320
PCRRXNS=3
POOLS=3
ARTHMINPROB=0.8
echo "Sumaclust similarity percentage is" ${SUMASIM}"%."
HOMEFOLDER="/Users/Negorashi2011/Xiaoyangmiseqdata/MiSeq_20170410/300Test/"  # do not have a ~ in the path
echo "Home folder is" ${HOMEFOLDER}
SEQS="data/seqs/"
ANALYSIS="analysis/"
DAME="/usr/local/bin/DAMe/bin/"

# install assert.sh
# cd ${HOMEFOLDER}scripts
# git clone https://github.com/torokmark/assert.sh.git
cd ${HOMEFOLDER}scripts/assert.sh/; source './assert.sh'  # assert.sh code, the foler for which i've downloaded into my scripts/ folder.  See https://github.com/torokmark/assert.sh for syntax

cd ${HOMEFOLDER}${SEQS} # cd into the sequence folder

#### Read in sample list and make a bash array of the sample names (e.g. A1_S1)
find * -maxdepth 0 -name "*_L001_R1_001.fastq.gz" > samplelist.txt  # find all files ending with _L001_R1_001_fastq.gz
sample_info=samplelist.txt # put samplelist.txt into variable
sample_names=($(cut -f 1 "${sample_info}" | uniq)) # convert variable to array this way
# echo "${sample_names[@]}" # to echo all array elements
# echo ${sample_names[0]} # to echo first array element
echo "There are" "${#sample_names[@]}" "samples that will be processed:  " "${sample_names[@]}" # echo number of elements in the array

# To run a loop interactively, select the entire loop and send to terminal.  Don't ctrl-Enter each line because this can send a command N times, as many lines as the command covers. So if the command wraps over 2 lines, ctrl-Entering each line sends the whole command twice, causing the command to be run twice per loop.

#### Use AdapterRemoval to trim Illumina sequencing adapters. (there are a few adapters still left in the raw data)
# brew install adapterremoval
# or download from https://github.com/MikkelSchubert/adapterremoval

# loop over all samples and remove any adapters
for sample in "${sample_names[@]}"  # ${sample_names[@]} is the full bash array.  So loop over all samples
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
for sample in "${sample_names[@]}"  # ${sample_names[@]} is the full bash array.  So loop over all samples
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
for sample in "${sample_names[@]}"  # ${sample_names[@]} is the full bash array.  So loop over all samples
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
for sample in "${sample_names[@]}"  # ${sample_names[@]} is the full bash array.  So loop over all samples
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
for sample in "${sample_names[@]}"  # ${sample_names[@]} is the full bash array.  So loop over all samples
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
for sample in "${sample_names[@]}"  # "${sample_names[@]}" is the full bash array.  So loop over all samples
do
    sample_prefix="$( basename ${sample} "_L001_R1_001.fastq.gz")"  # e.g. A1_S1 is a sample_prefix
    echo ${sample_prefix}
    mkdir folder_${sample_prefix}
    mv sickle_cor_panda_${sample_prefix}.fastq.gz folder_${sample_prefix}/
done


# 2. SORT  Sort through each fastq file and determine how many of each tag pair is in each fastq file. Each fastq file takes < 1 min
for sample in "${sample_names[@]}"  # ${sample_names[@]} is the full bash array.  So loop over all samples
do
              sample_prefix="$( basename $sample "_L001_R1_001.fastq.gz")"  # e.g. A1_S1 is a sample_prefix
              cd ${HOMEFOLDER}data/seqs # starting point
              echo ${sample_prefix}
              cd folder_${sample_prefix}
              python ${DAME}sort.py -fq sickle_cor_panda_${sample_prefix}.fastq.gz -p ${HOMEFOLDER}data/Primers_COILeray.txt -t ${HOMEFOLDER}data/Tags_300test_COI.txt
              # ls -lhS > sizesort_folder_${sample_prefix}.txt # quick check if the twin tag files are the largest. # sort by size.  The largest files should all be twin tag files (e.g. Tag1_Tag1.txt). Superseded by splitSummaryByPSInfo.py
done


# 3. Place PCR replicates (1, 2, 3) of the same experiment (e.g. A, B, ...) in the same folder and rename to pool{1,2,3}
# The three PCR replicate folders (on which sort.py was run) have to be in the same folder (e.g. 'folderA') and named 'pool1', 'pool2', and 'pool3'. No other pool folders

cd ${HOMEFOLDER}data/seqs

for sample in "${sample_names[@]}"  # ${sample_names[@]} is the full bash array.  So loop over all samples
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


# original loop without if then for mkdir
# for sample in "${sample_names[@]}"  # ${sample_names[@]} is the full bash array.  So loop over all samples
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

# echo "${sample_prefix}" | cut -c 1-2 # selects out the first and second character from sample_prefix
# echo "${sample_prefix}" | grep -o '^\D' # selects first character, using grep


# 3.0.1. code added 20171126:  replace Xiaoyangmiseqdata/MiSeq_20170410/300Test/data/seqs/FolderA/pool2/Tag27_Tag27.txt with Xiaoyangmiseqdata/MiSeq_20171121/data/seqs/folderA/pool2/Tag27_Tag27.txt
# the Tag27_Tag27.txt file in MiSeq_20170410/ is the re-run library, since Tag27_Tag27 (A2_mmmmbody) failed in the MiSeq_20170410 run

# cd ${HOMEFOLDER}data/seqs
# cp folderA folderA_original # make a copy in case something goes wrong
# cp ~/Xiaoyangmiseqdata/MiSeq_20171121/data/seqs/folderA/pool2/Tag27_Tag27.txt ${HOMEFOLDER}data/seqs/folderA/pool2/
## I also by-hand changed the line in MiSeq_20170410/300Test/data/seqs/FolderA/pool2/SummaryCounts.txt to reflect the updated sequence totals from the MiSeq_20171121 run.


# 3.1 Sort SummaryCounts.txt to SummaryCountsSorted.txt (Bohmann code)
cd ${HOMEFOLDER}data/seqs
for sample in "${sample_names[@]}"  # ${sample_names[@]} is the full bash array.  So loop over all samples
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
cd ${HOMEFOLDER}data/seqs
sample_libs=($(cat samplelist.txt | cut -c 1 | uniq))  # cut out all but the first letter of each filename and keep unique values, samplelist.txt should already exist
echo "There are" "${#sample_libs[@]}" "samples that will be processed:"  "${sample_libs[@]}" # echo number and name of elements in the array

# can replace the below with parallel syntax
# parallel echo "Moved SummaryCounts_sorted_Folder{1}_Pool{2}.txt" ::: ${sample_libs[@]} ::: `seq 1 ${POOLS}`

cd ${HOMEFOLDER}data/seqs
for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
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

for sample in "${sample_names[@]}"  # ${sample_names[@]} is the full bash array.  So loop over all samples
do
              sample_prefix="$( basename $sample "_L001_R1_001.fastq.gz")"
              cd ${HOMEFOLDER}data/seqs
              sample_prefix_prefix="$(echo "${sample_prefix}" | cut -c 1)"  # e.g. A is a sample_prefix_prefix
              sample_prefix_pool="$(echo "${sample_prefix}" | cut -c 2)"  # e.g. 1 is a sample_prefix_pool
              # cd folder${sample_prefix_prefix}/pool${sample_prefix_pool}  # cd into each pool (e.g. folderA/pool1)
              cd folder${sample_prefix_prefix} # cd into each folder (e.g. folderA)
              python ${DAME}splitSummaryByPSInfo.py -p ${HOMEFOLDER}data/PSinfo_300test_COI${sample_prefix_prefix}.txt -l ${sample_prefix_pool} -s pool${sample_prefix_pool}/SummaryCounts.txt -o splitSummaryByPSInfo_Folder${sample_prefix_prefix}_Pool${sample_prefix_pool}.txt
              date
done


# 3.3 Run scripts/heatmap.R on the different pools.

Rscript --vanilla --verbose ${HOMEFOLDER}scripts/heatmap.R

# then move heatmaps from inside pool folders into sample Folders
# should replace this with:  find | xargs mv
cd ${HOMEFOLDER}data/seqs
for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
              for pool in `seq 1 ${POOLS}`
              do
                            mv folder${sample}/pool${pool}/heatmap_Folder${sample}_Pool${pool}.pdf folder${sample}
              done
done


# 4  Filter

# 4.1 RSI test PCR replicate similarity.  First filter.py at -y 1 -t 1, to keep all sequences, then run RSI.py on the pools
              ### This step is slow (~ 1.3 hrs per library on my MacBookPro).  I run multiple jobs at once using parallel.

cd ${HOMEFOLDER}data/seqs
# python ${DAME}filter.py -h
# MINPCR_1=1 # min number of PCRs that a sequence has to appear in
# MINREADS_1=1 # min number of copies per sequence per PCR
# confirm MINPCR and MINREADS values
echo "For RSI analysis, all sequences are kept: appear in just ${MINPCR_1} PCR, with just ${MINREADS_1} read per PCR."

#### This bit of code to make sample_libs[] exists above too (SummaryCountsSorted.txt).  Running here again too, in case, i'm running this loop independently
#### Read in sample list and make a bash array of the sample libraries (A, B, C, D, E, F)
# find * -maxdepth 0 -name "*_L001_R1_001.fastq.gz" > samplelist.txt  # find all files ending with _L001_R1_001_fastq.gz
sample_libs=($(cat samplelist.txt | cut -c 1 | uniq))  # cut out all but the first letter of each filename and keep unique values, samplelist.txt should already exist
# echo ${sample_libs[1]} # to echo first array element
echo "There are" "${#sample_libs[@]}" "samples that will be processed:"  "${sample_libs[@]}" # echo number and name of elements in the array

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

# mkdir in parallel
parallel mkdir folder{1}/Filter_min${MINPCR_1}PCRs_min${MINREADS_1}copies_{1} ::: "${sample_libs[@]}"
rm filter1_1_commands.txt # ensure no filter1_1_commands.txt file is present
# create a list of commands with the correct arguments, using echo "command syntax"
for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
              echo "cd folder${sample}; \
              python ${DAME}filter.py -psInfo ${HOMEFOLDER}data/PSinfo_300test_COI${sample}.txt -x ${PCRRXNS} -y ${MINPCR_1} -p ${POOLS} -t ${MINREADS_1} -l ${MINLEN} -o Filter_min${MINPCR_1}PCRs_min${MINREADS_1}copies_${sample}; \
              python ${DAME}RSI.py --explicit -o RSI_output_${sample}.txt Filter_min${MINPCR_1}PCRs_min${MINREADS_1}copies_${sample}/Comparisons_${PCRRXNS}PCRs.txt" \
              >> filter1_1_commands.txt
done
# run parallel --dryrun to see the commands that will be run, without actually running them.
              parallel -k --dryrun :::: filter1_1_commands.txt
# run the command for real
              parallel --jobs 3 -k :::: filter1_1_commands.txt  # parallel :::: filter1_1_commands.txt means that the commands come from filter1_1_commands.txt
rm filter1_1_commands.txt # ensure no filter1_1_commands.txt file is present

# alternative parallel syntax, using composed commands (combined commands wrapped in single quotes).  still needs to be tested.
# parallel --dryrun 'cd ${HOMEFOLDER}data/seqs/folder{1}/pool{2}; ls Tag* | wc -l' ::: ${sample_libs[@]} ::: `seq 1 ${POOLS}` # example
# parallel --dryrun --jobs 3 -k 'cd folder{}; python ${DAME}filter.py -psInfo ${HOMEFOLDER}data/PSinfo_300test_COI{}.txt -x ${PCRRXNS} -y ${MINPCR_1} -p ${POOLS} -t ${MINREADS_1} -l ${MINLEN} -o Filter_min${MINPCR_1}PCRs_min${MINREADS_1}copies_{}; python ${DAME}RSI.py --explicit -o RSI_output_{}.txt Filter_min${MINPCR_1}PCRs_min${MINREADS_1}copies_{}/Comparisons_${PCRRXNS}PCRs.txt' ::: ${sample_libs[@]}


# non-parallel method.  analyses one folder at a time in a loop.
# cd ${HOMEFOLDER}data/seqs
# for sample in ${sample_libs[@]}  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
# do
#               cd ${HOMEFOLDER}data/seqs
#               cd folder${sample} # cd into folderA,B,C,D,E,F
#               mkdir Filter_min${MINPCR_1}PCRs_min${MINREADS_1}copies_${sample}
#               date
#
#               # filter.py
#               python ${DAME}filter.py -psInfo ${HOMEFOLDER}data/PSinfo_300test_COI${sample}.txt -x ${PCRRXNS} -y ${MINPCR_1} -p ${POOLS} -t ${MINREADS_1} -l ${MINLEN} -o Filter_min${MINPCR_1}PCRs_min${MINREADS_1}copies_${sample}
#
#               # RSI.py
#               python ${DAME}RSI.py --explicit -o RSI_output_${sample}.txt Filter_min${MINPCR_1}PCRs_min${MINREADS_1}copies_${sample}/Comparisons_${PCRRXNS}PCRs.txt
# done


# 4.2 plotLengthFreqMetrics_perSample.py and reads per sample per pool
# python ${DAME}plotLengthFreqMetrics_perSample.py -h
cd ${HOMEFOLDER}data/seqs

for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
              cd ${HOMEFOLDER}data/seqs
              cd folder${sample}/Filter_min${MINPCR_1}PCRs_min${MINREADS_1}copies_${sample} # cd into folderA,B,C,D,E,F/filter_output
              echo "Analysing: folder"${sample}/Filter_min${MINPCR_1}PCRs_min${MINREADS_1}copies_${sample}/
              # plotLengthFreqMetrics_perSample.py
              python ${DAME}plotLengthFreqMetrics_perSample.py -f FilteredReads.fna -p ${HOMEFOLDER}data/PSinfo_300test_COI${sample}.txt -n 3
done

# 4.2.1 change header lines to sumaclust format, remove chimeras using vsearch

# python ${DAME}convertToUSearch.py -h
# MINPCR=2 # these commands are to make it possible to process multiple filter.py outputs, which are saved in different Filter_min folders
# MINREADS=4

# installation of gsed (gsed == GNU version of sed == Linux version of sed)
# brew install gnu-sed

# vsearch uchime version, a few seconds per library when applied to FilteredReads.forsumaclust.fna
for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
              cd ${HOMEFOLDER}data/seqs/folder${sample}/Filter_min${MINPCR_1}PCRs_min${MINREADS_1}copies_${sample}
              python ${DAME}convertToUSearch.py -i FilteredReads.fna -lmin 300 -lmax 320
              gsed 's/ count/;size/' FilteredReads.forsumaclust.fna > FilteredReads.forvsearch.fna
              vsearch --sortbysize FilteredReads.forvsearch.fna --output FilteredReads.forvsearch_sorted.fna
              vsearch --uchime_denovo FilteredReads.forvsearch_sorted.fna --nonchimeras FilteredReads.forvsearch_sorted_nochimeras.fna
              gsed 's/;size/ count/' FilteredReads.forvsearch_sorted_nochimeras.fna > FilteredReads.forsumaclust.nochimeras.fna
done

# parallel vsearch uchime version, UNTESTED
# parallel --jobs 6 "cd ${HOMEFOLDER}data/seqs/folder{1}/Filter_min${MINPCR_1}PCRs_min${MINREADS_1}copies_{1}; python ${DAME}convertToUSearch.py -i FilteredReads.fna -lmin 300 -lmax 320;  gsed 's/ count/;size/' FilteredReads.forsumaclust.fna > FilteredReads.forvsearch.fna; vsearch --sortbysize FilteredReads.forvsearch.fna --output FilteredReads.forvsearch_sorted.fna; vsearch --uchime_denovo FilteredReads.forvsearch_sorted.fna --nonchimeras FilteredReads.forvsearch_sorted_nochimeras.fna; gsed 's/;size/ count/' FilteredReads.forvsearch_sorted_nochimeras.fna > FilteredReads.forsumaclust.nochimeras.fna" ::: "${sample_libs[@]}"

# remove vsearch uchime working files
for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
              cd ${HOMEFOLDER}data/seqs/folder${sample}/Filter_min${MINPCR_1}PCRs_min${MINREADS_1}copies_${sample}
              rm FilteredReads.forsumaclust.fna
              rm FilteredReads.forvsearch.fna
              rm FilteredReads.forvsearch_sorted.fna
              rm FilteredReads.forvsearch_sorted_nochimeras.fna
              rm Comparisons_${PCRRXNS}PCRs.fasta # large file that doesn't need to be kept
              rm Comparisons_${PCRRXNS}PCRs.txt # large file that doesn't need to be kept
done


# 4.2.2. Sumaclust clustering and convert Sumaclust output to table format

# sumaclust download from: https://git.metabarcoding.org/obitools/sumaclust/wikis/home
# wget https://git.metabarcoding.org/obitools/sumaclust/uploads/69f757c42f2cd45212c587e87c75a00f/sumaclust_v1.0.20.tar.gz
# tar -zxvf sumaclust_v1.0.20.tar.gz
# cd sumaclust_v1.0.20/
# make CC=clang # disables OpenMP, which isn't on macOS
# mv sumaclust /usr/local/bin/

# ~/src/sumaclust_v1.0.20/sumaclust -h
# python ${DAME}tabulateSumaclust.py -h

# 97% sumaclust
echo ${SUMASIM} # confirm the similarity value
SUMASIM=97 # if there is no SUMASIM value
echo ${SUMASIM} # confirm the similarity value
cd ${HOMEFOLDER}
for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
     cd ${HOMEFOLDER}data/seqs/folder${sample}/Filter_min${MINPCR_1}PCRs_min${MINREADS_1}copies_${sample}
     sumaclust -t .${SUMASIM} -e FilteredReads.forsumaclust.nochimeras.fna > OTUs_${SUMASIM}_sumaclust.fna
     python ${DAME}tabulateSumaclust.py -i OTUs_${SUMASIM}_sumaclust.fna -o table_300test_${sample}_${SUMASIM}.txt -blast
     mv table_300test_${sample}_${SUMASIM}.txt.blast.txt table_300test_${sample}_${SUMASIM}.fas # chg filename suffix, UNTESTED LINE
done


# create LULU matchfiles
for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
     cd ${HOMEFOLDER}data/seqs/folder${sample}/Filter_min${MINPCR_1}PCRs_min${MINREADS_1}copies_${sample}
     # remove the last column, which holds the
     vsearch --usearch_global table_300test_${sample}_${SUMASIM}.fas --db table_300test_${sample}_${SUMASIM}.fas --self --id .84 --iddef 1 --userout match_list_${sample}${SUMASIM}.txt -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10

     # makeblastdb -in table_300test_${sample}_${sim}_Arthropoda.fas -parse_seqids -dbtype nucl
     # blastn -db table_300test_${sample}_${sim}_Arthropoda.fas -outfmt '6 qseqid sseqid pident' -out match_list_${sample}${sim}.txt -qcov_hsp_perc 80 -perc_identity 84 -query table_300test_${sample}_${sim}_Arthropoda.fas
done


# 4.2.3. Move sumaclust results to ${HOMEFOLDER}/analysis
for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
              cd ${HOMEFOLDER}${ANALYSIS}
              # make folder to hold results
              if [ ! -d OTU_transient_results ] # if directory OTU_transient_results does not exist
              then
                  mkdir OTU_transient_results/
                  mkdir OTU_transient_results/OTU_tables/
              fi
              mv ${HOMEFOLDER}${SEQS}folder${sample}/Filter_min${MINPCR_1}PCRs_min${MINREADS_1}copies_${sample} ${HOMEFOLDER}${ANALYSIS}OTU_transient_results/Filter_min${MINPCR_1}PCRs_min${MINREADS_1}copies_${sample}
              cd ${HOMEFOLDER}${ANALYSIS}OTU_transient_results/Filter_min${MINPCR_1}PCRs_min${MINREADS_1}copies_${sample}
              cp table_300test_${sample}_*.txt ${HOMEFOLDER}${ANALYSIS}OTU_transient_results/OTU_tables/ # copy OTU tables to OTU_tables folder
              cp table_300test_${sample}_*.fas ${HOMEFOLDER}${ANALYSIS}OTU_transient_results/OTU_tables/ # copy OTU seqs to OTU_tables folder
              cp match_list_${sample}${SUMASIM}.txt ${HOMEFOLDER}${ANALYSIS}OTU_transient_results/OTU_tables/ # copy matchlist to OTU_tables folder
done

mv ${HOMEFOLDER}${ANALYSIS}OTU_transient_results/ ${HOMEFOLDER}${ANALYSIS}OTUs_min${MINPCR_1}PCRs_min${MINREADS_1}copies_"$(date +%F_time-%H%M)"/


# 4.2.4.  Run LULU code here to collapse the OTU table

OTUTABLEFOLDER="OTUs_min1PCRs_min1copies_2017-11-30_time-2148" # this needs to be set by hand because the time stamp is unpredictable

Rscript --vanilla --verbose ${HOMEFOLDER}scripts/LULU.R ${OTUTABLEFOLDER}
rm ${HOMEFOLDER}${ANALYSIS}/${OTUTABLEFOLDER}/OTU_tables/lulu.log_*

# filter fasta files to keep only the OTUs that are kept by LULU
for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
     for sim in `seq 97 97`  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
     do
          cd ${HOMEFOLDER}${ANALYSIS}/${OTUTABLEFOLDER}/OTU_tables
          seqtk subseq table_300test_${sample}_${sim}.fas <(cut -f 1 table_300test_${sample}_${sim}_lulu.txt) > table_300test_${sample}_${sim}_lulu.fas
               # <(cut -f 1 table_300test_${sample}_${sim}_lulu.txt) produces the list of sequences to keep
     done
done

# Upload to hpc to do taxonomic assignment with midori

# Download the RDP output files back to the OTU_tables folder

# 4.2.5. Filter out non-Arthropoda from RDP assignment table.  Keep only Arthropoda with prob >= ARTHMINPROB (set to 0.80).
# Filter out non-Arthropoda OTUs from OTU representative sequences fasta file
for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
     for sim in `seq 97 97`  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
     do
          cd ${HOMEFOLDER}${ANALYSIS}/${OTUTABLEFOLDER}/OTU_tables
          awk -v arthmin=${ARTHMINPROB} '$8 ~ /Arthropoda/ && $10 >= arthmin { print }' table_300test_${sample}_${sim}.RDPmidori_lulu.txt > table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt
          gsed -E 's/\t\t/\t/' table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt > table_300test_${sample}_${sim}.RDPmidori_Arthropoda_nodbltab.txt
          mv table_300test_${sample}_${sim}.RDPmidori_Arthropoda_nodbltab.txt table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt
          seqtk subseq table_300test_${sample}_${sim}.fas <(cut -f 1 table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt) > table_300test_${sample}_${sim}_Arthropoda.fas
               # <(cut -f 1 table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt) produces the list of sequences to keep
     done
done

# 4.2.6. checking that the right number of OTUs has been removed from each file
for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
     for sim in `seq 97 97`  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
     do
          echo "OTU tables"
          wc -l table_300test_${sample}_${sim}.RDPmidori_lulu.txt
          wc -l table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt
          echo "fasta files"
          grep ">" table_300test_${sample}_${sim}_lulu.fas | wc -l
          grep ">" table_300test_${sample}_${sim}_Arthropoda.fas | wc -l
     done
done

# Hemiptera;Caliscelidae;Bruchomorpha is missing its family in the Midori database. So the identifcation omits the family name. This prevents R from inputting the OTU table. Eventually, I need to change the MIDORI database and retrain it.
# grep "Bruchomorpha" MIDORI_UNIQUE_1.1.1_COI_RDP.fasta
# sed 's/Hemiptera;Bruchomorpha/Hemiptera;Caliscelidae;Bruchomorpha/' MIDORI_UNIQUE_1.1_COI_RDP.fasta > MIDORI_UNIQUE_1.1.1_COI_RDP.fasta

# DANGEROUS CODE:  RUN ONLY ONCE AFTER GENERATING THE RDP ARTHROPODA-ONLY TABLES, BECAUSE IF RUN MORE THAN ONCE, WILL INSERT THE NEW TAXONOMIC RANK (e.g. Caliscelidae family 0.5)  MORE THAN ONCE
# OTUTABLEFOLDER="OTUs_min2PCRs_min4copies_2017-11-26_time-1727"
# for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
# do
#      for sim in `seq 97 97`  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
#      do
#           cd ${HOMEFOLDER}${ANALYSIS}/${OTUTABLEFOLDER}/OTU_tables
#           gsed -E 's/Bruchomorpha\tgenus/Caliscelidae\tfamily\t0.50\t\Bruchomorpha\tgenus/' table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt > table_300test_${sample}_${sim}.RDPmidori_Arthropoda_family.txt
#           mv table_300test_${sample}_${sim}.RDPmidori_Arthropoda_family.txt table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt
#           echo "table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt"
#           grep "Bruchomorpha" table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt
#      done
# done

# DANGEROUS CODE:  RUN ONLY ONCE AFTER GENERATING THE RDP ARTHROPODA-ONLY TABLES, BECAUSE IF RUN MORE THAN ONCE, WILL INSERT THE NEW TAXONOMIC RANK (e.g. Caliscelidae family 0.5)  MORE THAN ONCE
# for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
# do
#      for sim in `seq 97 97`
#      do
#           cd ${HOMEFOLDER}${ANALYSIS}/${OTUTABLEFOLDER}/OTU_tables
#           gsed -E 's/Entomobryidae\tfamily/Entomobryomorpha\torder\t0.50\tEntomobryidae\tfamily/' table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt > table_300test_${sample}_${sim}.RDPmidori_Arthropoda_order.txt
#           mv table_300test_${sample}_${sim}.RDPmidori_Arthropoda_order.txt table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt
#           echo "table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt"
#           grep "Entomobryomorpha" table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt
#      done
# done

##### BLAST the final OTUs against the MTB reference sequences. Note that the MTB reference seqs are already clustered at 97% via sumaclust. There are 254 sumaclust 97% OTUs.

cd ${HOMEFOLDER}/data/MTB
# run once
# makeblastdb -in MTB_AllInputRefSeqs_20170726.fasta -dbtype nucl # to make the MTB ref dataset BLASTABLE

echo $SUMASIM # should be 97

# take the top hit from blast output. visual inspection of the blast hits with similarity < 0.98 shows that the low-similarity hits are echo OTUs (echo OTUs are small OTUs that are similar to large OTUs that hit an MTB sequence at ~100% similarity)
parallel --jobs 4 "blastn -db ${HOMEFOLDER}data/MTB/MTB_AllInputRefSeqs_20170726.fasta -query ${HOMEFOLDER}analysis/${OTUTABLEFOLDER}/OTU_tables/table_300test_{1}_${SUMASIM}_Arthropoda.fas -num_threads 2 -evalue 1e-10 -max_target_seqs 1 -outfmt 6 -out ${HOMEFOLDER}analysis/${OTUTABLEFOLDER}/OTU_tables/table_300test_{1}_${SUMASIM}_Arthropoda.blastnMTB.txt" ::: A B C D E F

# single library code
# EXPERIMENT=F
# blastn -db ${HOMEFOLDER}/data/MTB/MTB_AllInputRefSeqs_20170726.fasta -query ${HOMEFOLDER}/analysis/${OTUTABLEFOLDER}/OTU_tables/table_300test_${EXPERIMENT}_${SUMASIM}_Arthropoda.fas -num_threads 6 -evalue 1e-10 -max_target_seqs 1 -outfmt 6 -out ${HOMEFOLDER}/analysis/${OTUTABLEFOLDER}/OTU_tables/table_300test_${EXPERIMENT}_${SUMASIM}_Arthropoda.blastnMTB.txt

# outfmt 6 column headings
# qseqid	sseqid	pident	length	mismatch	gapopen	qstart	qend	sstart	send	evalue	bitscore

# now run Dropout analyses in R


####################################################################################################
## After consideration of the negative controls, summary counts, and the heatmaps, choose thresholds for filtering
# Will probably have to run a few times (+ OTU clustering) with different values of minPCR and minReads to find values that result in no OTUs from the extraction blank and PCR blank negative controls
####################################################################################################

# 4.3 Filter - Filtering reads with minPCR and minReads/PCR thresholds, min 300 bp length.
              ### This step is slow (~ 1.3 hr per library, but i run up to 4 jobs in parallel on my macbookpro).

#### If I want to re-run filter.py with different thresholds, I set new values of MINPCR & MINREADS and start from here.
MINPCR=2 # min number of PCRs that a sequence has to appear in
MINREADS=4 # min number of copies per sequence per PCR
# confirm MINPCR and MINREADS values
echo "Each (unique) sequence must appear in at least ${MINPCR} PCRs, with at least ${MINREADS} reads per PCR."

cd ${HOMEFOLDER}data/seqs
# python ${DAME}filter.py -h

#### Read in sample list and make a bash array of the sample libraries (A, B, C, D, E, F)
# find * -maxdepth 0 -name "*_L001_R1_001.fastq.gz" > samplelist.txt  # find all files ending with _L001_R1_001_fastq.gz
sample_libs=($(cat samplelist.txt | cut -c 1 | uniq))  # cut out all but the first letter of each filename and keep unique values, samplelist.txt should already exist
# echo ${sample_libs[1]} # to echo first array element
echo "There are" "${#sample_libs[@]}" "samples that will be processed:" "${sample_libs[@]}" # echo number and name of elements in the array

cd ${HOMEFOLDER}data/seqs
# mkdir in parallel
parallel mkdir folder{1}/Filter_min${MINPCR}PCRs_min${MINREADS}copies_{1} ::: "${sample_libs[@]}"
rm filter_commands.txt # ensure no filter_commands.txt file is present
# create a list of commands with the correct arguments
for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
              echo "cd folder${sample}; \
              python ${DAME}filter.py -psInfo ${HOMEFOLDER}data/PSinfo_300test_COI${sample}.txt -x ${PCRRXNS} -y ${MINPCR} -p ${POOLS} -t ${MINREADS} -l ${MINLEN} -o Filter_min${MINPCR}PCRs_min${MINREADS}copies_${sample}" \
              >> filter_commands.txt
done
# run parallel --dryrun to see the commands that will be run, without actually running them.
parallel --jobs 4 --eta -k :::: filter_commands.txt  # parallel :::: filter_commands.txt means that the commands come from filter_commands.txt.  --eta estimates how long till the jobs finish, using individual job times

rm filter_commands.txt # ensure no filter_commands.txt file is present

# alternative parallel syntax, using composed commands.  still needs to be tested.
# parallel --dryrun --jobs 3 -k 'cd folder{}; python ${DAME}filter.py -psInfo ${HOMEFOLDER}data/PSinfo_300test_COI{}.txt -x ${PCRRXNS} -y ${MINPCR} -p ${POOLS} -t ${MINREADS} -l ${MINLEN} -o Filter_min${MINPCR}PCRs_min${MINREADS}copies_{}' ::: ${sample_libs[@]}

# alternative loop method.  analyses one folder at a time
              # for sample in ${sample_libs[@]}  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
              # do
              #               cd ${HOMEFOLDER}data/seqs
              #               cd folder${sample} # cd into folderA,B,C,D,E,F
              #               mkdir Filter_min${MINPCR}PCRs_min${MINREADS}copies_${sample}
              #               date
              #               python ${DAME}filter.py -psInfo ${HOMEFOLDER}data/PSinfo_300test_COI${sample}.txt -x ${PCRRXNS} -y ${MINPCR} -p ${POOLS} -t ${MINREADS} -l ${MINLEN} -o Filter_min${MINPCR}PCRs_min${MINREADS}copies_${sample}
              # done

# python ${DAME}filter.py -h
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
# python ${DAME}plotLengthFreqMetrics_perSample.py -h

cd ${HOMEFOLDER}data/seqs

for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
              cd ${HOMEFOLDER}data/seqs
              cd folder${sample}/Filter_min${MINPCR}PCRs_min${MINREADS}copies_${sample} # cd into folderA,B,C,D,E,F/filter_output
              echo "Analysing: folder"${sample}/Filter_min${MINPCR}PCRs_min${MINREADS}copies_${sample}/
              # plotLengthFreqMetrics_perSample.py
              python ${DAME}plotLengthFreqMetrics_perSample.py -f FilteredReads.fna -p ${HOMEFOLDER}data/PSinfo_300test_COI${sample}.txt -n 3
done


# 5. Prepare the FilteredReads.fna file for Sumaclust clustering: analyse pairwise similarities, change the header lines, and remove chimeras.

# install sumatra
              # cd ~/src
              # wget https://git.metabarcoding.org/obitools/sumatra/uploads/251020bbbd6c6595cb9fce6077e29952/sumatra_v1.0.20.tar.gz
              # tar -zxvf sumatra_v1.0.20.tar.gz
              # cd sumatra_v1.0.20/
              # make CC=clang # disables OpenMP, which isn't on macOS
              # mv sumatra /usr/local/bin

# 5.1 assessing clustering parameters
# This takes quite a long time on large fna files, only feasible on FilteredReads.fna that have been filtered to a few MB.  Takes about 5 mins per fna file on my macbookpro, after filtering down to min2PCRs_min4copies
cd ${HOMEFOLDER}data/seqs/
# MINPCR=2 # these commands are to make it possible to process multiple filter.py outputs, which are saved in different Filter_min folders
# MINREADS=4
echo "Analysing filter.py output where each (unique) sequence appeared in ≥ ${MINPCR} PCRs, with ≥ ${MINREADS} reads per PCR."

# for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
# do
#               cd ${HOMEFOLDER}data/seqs/folder${sample}/Filter_min${MINPCR}PCRs_min${MINREADS}copies_${sample}
#               python ${DAME}assessClusteringParameters.py -i FilteredReads.fna -mint 0.8 -minR 0.6 -step 0.05 -t 4 -o COIclusterassess_mint08_minR06_step005_Filter_min${MINPCR}PCRs_min${MINREADS}copies_${sample}.pdf
# done

# untested code using gnu parallel, but it should work.  needs to be tested.
rm cluster_commands.txt # ensure no cluster_commands.txt file is present
# create a list of commands with the correct arguments
for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
              echo "cd ${HOMEFOLDER}data/seqs/folder${sample}/Filter_min${MINPCR}PCRs_min${MINREADS}copies_${sample}; \
              python ${DAME}assessClusteringParameters.py -i FilteredReads.fna -mint 0.8 -minR 0.6 -step 0.05 -t 4 -o COIclusterassess_mint08_minR06_step005_Filter_min${MINPCR}PCRs_min${MINREADS}copies_${sample}.pdf" \
              >> cluster_commands.txt
done
# run parallel --dryrun to see the commands that will be run, without actually running them.
parallel --jobs 4 --eta -k :::: cluster_commands.txt  # parallel :::: cluster_commands.txt means that the commands come from cluster_commands.txt.  --eta estimates how long till the jobs finish, using individual job times

rm cluster_commands.txt # ensure no cluster_commands.txt file is present

              # python ${DAME}assessClusteringParameters.py -h
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

# 5.2 change header lines to sumaclust format, remove chimeras using vsearch

# python ${DAME}convertToUSearch.py -h
# MINPCR=2 # these commands are to make it possible to process multiple filter.py outputs, which are saved in different Filter_min folders
# MINREADS=4

# installation of gsed (gsed == GNU version of sed == Linux version of sed)
# brew install gnu-sed

# vsearch uchime version, a few seconds per library when applied to FilteredReads.forsumaclust.fna
for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
              cd ${HOMEFOLDER}data/seqs/folder${sample}/Filter_min${MINPCR}PCRs_min${MINREADS}copies_${sample}
              python ${DAME}convertToUSearch.py -i FilteredReads.fna -lmin 300 -lmax 320
              gsed 's/ count/;size/' FilteredReads.forsumaclust.fna > FilteredReads.forvsearch.fna
              vsearch --sortbysize FilteredReads.forvsearch.fna --output FilteredReads.forvsearch_sorted.fna
              vsearch --uchime_denovo FilteredReads.forvsearch_sorted.fna --nonchimeras FilteredReads.forvsearch_sorted_nochimeras.fna
              gsed 's/;size/ count/' FilteredReads.forvsearch_sorted_nochimeras.fna > FilteredReads.forsumaclust.nochimeras.fna
done

# remove vsearch uchime working files
for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
              cd ${HOMEFOLDER}data/seqs/folder${sample}/Filter_min${MINPCR}PCRs_min${MINREADS}copies_${sample}
              rm FilteredReads.forsumaclust.fna
              rm FilteredReads.forvsearch.fna
              rm FilteredReads.forvsearch_sorted.fna
              rm FilteredReads.forvsearch_sorted_nochimeras.fna
done


# Typically will find < 25 unique reads that are chimeras (<0.1%), but some of these reads would turn into their own OTUs

# usearch uchime version
# for sample in ${sample_libs[@]}  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
# do
#               cd ${HOMEFOLDER}data/seqs/folder${sample}/Filter_min${MINPCR}PCRs_min${MINREADS}copies_${sample}
#               python ${DAME}convertToUSearch.py -i FilteredReads.fna -lmin 300 -lmax 320
#               gsed 's/ count/;size/' FilteredReads.forsumaclust.fna > FilteredReads.foruchime.fna
#               usearch9 -sortbysize FilteredReads.foruchime.fna -fastaout FilteredReads.foruchime_sorted.fna
#               usearch9 -uchime2_denovo FilteredReads.foruchime_sorted.fna -nonchimeras FilteredReads.foruchime_sorted_nochimeras.fna
#               gsed 's/;size/ count/' FilteredReads.foruchime_sorted_nochimeras.fna > FilteredReads.forsumaclust.nochimeras.fna
#               # python ${DAME}convertToUSearch.py -i FilteredReads.fna -lmin 300 -lmax 320 -u
# done
# full usearch -uchime2_denovo command syntax
# usearch9 -uchime2_denovo FilteredReads.foruchime_sorted.fna -uchimeout out.txt -chimeras ch.fa -nonchimeras nonch.fa
# usearch10 -uchime3_denovo does not work

# 6. Sumaclust clustering and convert Sumaclust output to table format

# sumaclust download from: https://git.metabarcoding.org/obitools/sumaclust/wikis/home
# wget https://git.metabarcoding.org/obitools/sumaclust/uploads/69f757c42f2cd45212c587e87c75a00f/sumaclust_v1.0.20.tar.gz
# tar -zxvf sumaclust_v1.0.20.tar.gz
# cd sumaclust_v1.0.20/
# make CC=clang # disables OpenMP, which isn't on macOS
# mv sumaclust /usr/local/bin/

# ~/src/sumaclust_v1.0.20/sumaclust -h
# python ${DAME}tabulateSumaclust.py -h

# 96% sumaclust
echo ${SUMASIM} # confirm that there is a similarity value chosen
SUMASIM=96 # if there is no SUMASIM value
echo ${SUMASIM} # confirm that there is a similarity value chosen
cd ${HOMEFOLDER}
for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
              cd ${HOMEFOLDER}data/seqs/folder${sample}/Filter_min${MINPCR}PCRs_min${MINREADS}copies_${sample}
              sumaclust -t .${SUMASIM} -e FilteredReads.forsumaclust.nochimeras.fna > OTUs_${SUMASIM}_sumaclust.fna
              python ${DAME}tabulateSumaclust.py -i OTUs_${SUMASIM}_sumaclust.fna -o table_300test_${sample}_${SUMASIM}.txt -blast
              mv table_300test_${sample}_${SUMASIM}.txt.blast.txt table_300test_${sample}_${SUMASIM}.fas # chg filename suffix, UNTESTED LINE
done


# 97% sumaclust
# MINPCR=2 # these commands are to make it possible to prepare multiple filter.py outputs
# MINREADS=4
echo ${SUMASIM} # confirm the similarity value
SUMASIM=97 # if there is no SUMASIM value
echo ${SUMASIM} # confirm the similarity value
cd ${HOMEFOLDER}
for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
              cd ${HOMEFOLDER}data/seqs/folder${sample}/Filter_min${MINPCR}PCRs_min${MINREADS}copies_${sample}
              sumaclust -t .${SUMASIM} -e FilteredReads.forsumaclust.nochimeras.fna > OTUs_${SUMASIM}_sumaclust.fna
              python ${DAME}tabulateSumaclust.py -i OTUs_${SUMASIM}_sumaclust.fna -o table_300test_${sample}_${SUMASIM}.txt -blast
              mv table_300test_${sample}_${SUMASIM}.txt.blast.txt table_300test_${sample}_${SUMASIM}.fas # chg filename suffix, UNTESTED LINE
done

#### ***** Using CROP makes fewer OTUs, compared to sumaclust
# cd folderA/Filter_min2PCRs_min2copies_A/
# ~/src/CROP/CROP -i FilteredReads.forsumaclust.fna -o OTUs_CROP_97.out -s -e 2000 -b 185 -z 300 -r 0
              # Parms are:
              # -r = "This parameter controls the size of the clusters to be considered as “rare”.  –r 0 will yield the best accuracy
              # -z < 150,000/avg_seq_length = 150,000/320 bp = 469. let zl < 150000, so set z=300;
              # -b = block size = number_of_sequences/50 = 9227/50 = 185;  300Test/folderA/Filter_min2PCRs_min2copies_A/FilteredReads.forsumaclust.fna
              # -e = 10*block size = 1850.  2000 (default)
              # -s # 97% similarity
#  HOWEVER, CROP takes 5 minutes, and sumaclust takes 10 secs. Also, CROP now needs a script to convert the CROP OTU map (OTUs_CROP_97.out.cluster.list) and the OTU representative seq list (OTUs_CROP_97.out.cluster.fasta) into an OTU table that preserves the number of original Illumina reads per OTU-sample combination (e.g. table_300test_A_97.txt, which is post sumaclust + tabulateSumaclust.py + OTUs_97_sumaclust.fna)

# 7. Move sumaclust results to ${HOMEFOLDER}/analysis

# MINPCR=2; MINREADS=4 # if not running as bash

for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
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
              cp table_300test_${sample}_*.txt ${HOMEFOLDER}${ANALYSIS}OTU_transient_results/OTU_tables/ # copy OTU tables to OTU_tables folder
              cp table_300test_${sample}_*.fas ${HOMEFOLDER}${ANALYSIS}OTU_transient_results/OTU_tables/ # copy OTU seqs to OTU_tables folder, UNTESTED LINE
done

# change name of OTU_transient_results folder to include filter.py thresholds and timestamp
mv ${HOMEFOLDER}${ANALYSIS}OTU_transient_results/ ${HOMEFOLDER}${ANALYSIS}OTUs_min${MINPCR}PCRs_min${MINREADS}copies_"$(date +%F_time-%H%M)"/

# 8.0 create LULU matchfiles
OTUTABLEFOLDER="OTUs_min2PCRs_min4copies_2017-11-26_time-1727" # this needs to be set by hand because the time stamp is unpredictable

for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
     # cd ${HOMEFOLDER}data/seqs/folder${sample}/Filter_min${MINPCR}PCRs_min${MINREADS}copies_${sample}
     cd ${HOMEFOLDER}${ANALYSIS}/${OTUTABLEFOLDER}/OTU_tables
     # remove the last column, which holds the
     vsearch --usearch_global table_300test_${sample}_${SUMASIM}.fas --db table_300test_${sample}_${SUMASIM}.fas --self --id .84 --iddef 1 --userout match_list_${sample}${SUMASIM}.txt -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10

     # makeblastdb -in table_300test_${sample}_${sim}_Arthropoda.fas -parse_seqids -dbtype nucl
     # blastn -db table_300test_${sample}_${sim}_Arthropoda.fas -outfmt '6 qseqid sseqid pident' -out match_list_${sample}${sim}.txt -qcov_hsp_perc 80 -perc_identity 84 -query table_300test_${sample}_${sim}_Arthropoda.fas
done

# 8.1  Run LULU code here to collapse the OTU table

Rscript --vanilla --verbose ${HOMEFOLDER}scripts/LULU.R ${OTUTABLEFOLDER}
rm ${HOMEFOLDER}${ANALYSIS}/${OTUTABLEFOLDER}/OTU_tables/lulu.log_*

# filter fasta files to keep only the OTUs that are kept by LULU
for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
     for sim in `seq 97 97`  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
     do
          cd ${HOMEFOLDER}${ANALYSIS}/${OTUTABLEFOLDER}/OTU_tables
          seqtk subseq table_300test_${sample}_${sim}.fas <(cut -f 1 table_300test_${sample}_${sim}_lulu.txt) > table_300test_${sample}_${sim}_lulu.fas
               # <(cut -f 1 table_300test_${sample}_${sim}_lulu.txt) produces the list of sequences to keep
     done
done


# 8. After sumaclust, I upload the OTU fasta files to hpc.uea.ac.uk and assign taxonomies via RDP Classifier on the Midori database, and I filter out all non-Arthropoda (all in ~/midori/).  This has to run on the hpc.uea.ac.uk server because it needs ~18GB of RAM. Need to experiment on macOS, but otherwise can run remotely on hpc.

# Download the RDP output files back to the OTU_tables folder

# 9.0 Filter out non-Arthropoda from RDP assignment table.  Keep only Arthropoda with prob >= ARTHMINPROB (set to 0.80).
# Filter out non-Arthropoda OTUs from OTU representative sequences fasta file
for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
     for sim in `seq 97 97`  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
     do
          cd ${HOMEFOLDER}${ANALYSIS}/${OTUTABLEFOLDER}/OTU_tables
          awk -v arthmin=${ARTHMINPROB} '$8 ~ /Arthropoda/ && $10 >= arthmin { print }' table_300test_${sample}_${sim}.RDPmidori_lulu.txt > table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt
          gsed -E 's/\t\t/\t/' table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt > table_300test_${sample}_${sim}.RDPmidori_Arthropoda_nodbltab.txt
          mv table_300test_${sample}_${sim}.RDPmidori_Arthropoda_nodbltab.txt table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt
          seqtk subseq table_300test_${sample}_${sim}.fas <(cut -f 1 table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt) > table_300test_${sample}_${sim}_Arthropoda.fas
               # <(cut -f 1 table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt) produces the list of sequences to keep
     done
done

# 9.1. checking that the right number of OTUs has been removed from each file
for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
     for sim in `seq 97 97`  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
     do
          echo "OTU tables"
          wc -l table_300test_${sample}_${sim}.RDPmidori_lulu.txt
          wc -l table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt
          echo "fasta files"
          grep ">" table_300test_${sample}_${sim}_lulu.fas | wc -l
          grep ">" table_300test_${sample}_${sim}_Arthropoda.fas | wc -l
     done
done



# Hemiptera;Caliscelidae;Bruchomorpha is missing its family in the Midori database. So the identifcation omits the family name. This prevents R from inputting the OTU table. Eventually, I need to change the MIDORI database and retrain it.
# grep "Bruchomorpha" MIDORI_UNIQUE_1.1.1_COI_RDP.fasta
# sed 's/Hemiptera;Bruchomorpha/Hemiptera;Caliscelidae;Bruchomorpha/' MIDORI_UNIQUE_1.1_COI_RDP.fasta > MIDORI_UNIQUE_1.1.1_COI_RDP.fasta

# DANGEROUS CODE:  RUN ONLY ONCE AFTER GENERATING THE RDP ARTHROPODA-ONLY TABLES, BECAUSE IF RUN MORE THAN ONCE, WILL INSERT THE NEW TAXONOMIC RANK (e.g. Caliscelidae family 0.5)  MORE THAN ONCE
# OTUTABLEFOLDER="OTUs_min2PCRs_min4copies_2017-11-26_time-1727"
# for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
# do
#      for sim in `seq 97 97`  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
#      do
#           cd ${HOMEFOLDER}${ANALYSIS}/${OTUTABLEFOLDER}/OTU_tables
#           gsed -E 's/Bruchomorpha\tgenus/Caliscelidae\tfamily\t0.50\t\Bruchomorpha\tgenus/' table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt > table_300test_${sample}_${sim}.RDPmidori_Arthropoda_family.txt
#           mv table_300test_${sample}_${sim}.RDPmidori_Arthropoda_family.txt table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt
#           echo "table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt"
#           grep "Bruchomorpha" table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt
#      done
# done

# DANGEROUS CODE:  RUN ONLY ONCE AFTER GENERATING THE RDP ARTHROPODA-ONLY TABLES, BECAUSE IF RUN MORE THAN ONCE, WILL INSERT THE NEW TAXONOMIC RANK (e.g. Caliscelidae family 0.5)  MORE THAN ONCE
# for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
# do
#      for sim in `seq 97 97`
#      do
#           cd ${HOMEFOLDER}${ANALYSIS}/${OTUTABLEFOLDER}/OTU_tables
#           gsed -E 's/Hypsidia\tgenus/Drepanidae\tfamily\t0.50\tHypsidia\tgenus/' table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt > table_300test_${sample}_${sim}.RDPmidori_Arthropoda_order.txt
#           mv table_300test_${sample}_${sim}.RDPmidori_Arthropoda_order.txt table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt
#           echo "table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt"
#           grep "Drepanidae" table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt
#      done
# done

# DANGEROUS CODE:  RUN ONLY ONCE AFTER GENERATING THE RDP ARTHROPODA-ONLY TABLES, BECAUSE IF RUN MORE THAN ONCE, WILL INSERT THE NEW TAXONOMIC RANK (e.g. Caliscelidae family 0.5)  MORE THAN ONCE
# for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
# do
#      for sim in `seq 97 97`
#      do
#           cd ${HOMEFOLDER}${ANALYSIS}/${OTUTABLEFOLDER}/OTU_tables
#           gsed -E 's/Choroterpes\tgenus/Leptophlebiidae\tfamily\t0.50\tChoroterpes\tgenus/' table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt > table_300test_${sample}_${sim}.RDPmidori_Arthropoda_order.txt
#           mv table_300test_${sample}_${sim}.RDPmidori_Arthropoda_order.txt table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt
#           echo "table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt"
#           grep "Leptophlebiidae" table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt
#      done
# done

# DANGEROUS CODE:  RUN ONLY ONCE AFTER GENERATING THE RDP ARTHROPODA-ONLY TABLES, BECAUSE IF RUN MORE THAN ONCE, WILL INSERT THE NEW TAXONOMIC RANK (e.g. Caliscelidae family 0.5)  MORE THAN ONCE
# for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
# do
#      for sim in `seq 97 97`
#      do
#           cd ${HOMEFOLDER}${ANALYSIS}/${OTUTABLEFOLDER}/OTU_tables
#           gsed -E 's/Plectrotarsus\tgenus/Plectrotarsidae\tfamily\t0.50\tPlectrotarsus\tgenus/' table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt > table_300test_${sample}_${sim}.RDPmidori_Arthropoda_order.txt
#           mv table_300test_${sample}_${sim}.RDPmidori_Arthropoda_order.txt table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt
#           echo "table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt"
#           grep "Plectrotarsidae" table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt
#      done
# done

# DANGEROUS CODE:  RUN ONLY ONCE AFTER GENERATING THE RDP ARTHROPODA-ONLY TABLES, BECAUSE IF RUN MORE THAN ONCE, WILL INSERT THE NEW TAXONOMIC RANK (e.g. Caliscelidae family 0.5)  MORE THAN ONCE
# for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
# do
#      for sim in `seq 97 97`
#      do
#           cd ${HOMEFOLDER}${ANALYSIS}/${OTUTABLEFOLDER}/OTU_tables
#           gsed -E 's/Lepidotrichidae\tfamily/Zygentoma\torder\t0.50\tLepidotrichidae\tfamily/' table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt > table_300test_${sample}_${sim}.RDPmidori_Arthropoda_order.txt
#           mv table_300test_${sample}_${sim}.RDPmidori_Arthropoda_order.txt table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt
#           echo "table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt"
#           grep "Lepidotrichidae" table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt
#      done
# done

# DANGEROUS CODE:  RUN ONLY ONCE AFTER GENERATING THE RDP ARTHROPODA-ONLY TABLES, BECAUSE IF RUN MORE THAN ONCE, WILL INSERT THE NEW TAXONOMIC RANK (e.g. Caliscelidae family 0.5)  MORE THAN ONCE
# for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
# do
#      for sim in `seq 97 97`
#      do
#           cd ${HOMEFOLDER}${ANALYSIS}/${OTUTABLEFOLDER}/OTU_tables
#           gsed -E 's/Entomobryidae\tfamily/Entomobryomorpha\torder\t0.50\tEntomobryidae\tfamily/' table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt > table_300test_${sample}_${sim}.RDPmidori_Arthropoda_order.txt
#           mv table_300test_${sample}_${sim}.RDPmidori_Arthropoda_order.txt table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt
#           echo "table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt"
#           grep "Entomobryomorpha" table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt
#      done
# done

# DANGEROUS CODE:  RUN ONLY ONCE AFTER GENERATING THE RDP ARTHROPODA-ONLY TABLES, BECAUSE IF RUN MORE THAN ONCE, WILL INSERT THE NEW TAXONOMIC RANK (e.g. Caliscelidae family 0.5)  MORE THAN ONCE
# for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
# do
#      for sim in `seq 97 97`
#      do
#           cd ${HOMEFOLDER}${ANALYSIS}/${OTUTABLEFOLDER}/OTU_tables
#           gsed -E 's/Habrosyne\tgenus/Drepanidae\tfamily\t0.50\tHabrosyne\tgenus/' table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt > table_300test_${sample}_${sim}.RDPmidori_Arthropoda_order.txt
#           mv table_300test_${sample}_${sim}.RDPmidori_Arthropoda_order.txt table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt
#           echo "table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt"
#           grep "Drepanidae" table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt
#      done
# done

# DANGEROUS CODE:  RUN ONLY ONCE AFTER GENERATING THE RDP ARTHROPODA-ONLY TABLES, BECAUSE IF RUN MORE THAN ONCE, WILL INSERT THE NEW TAXONOMIC RANK (e.g. Caliscelidae family 0.5)  MORE THAN ONCE
# for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
# do
#      for sim in `seq 97 97`
#      do
#           cd ${HOMEFOLDER}${ANALYSIS}/${OTUTABLEFOLDER}
#           gsed -E 's/Coronididae\tfamily/Stomatopoda\torder\t0.50\tCoronididae\tfamily/' table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt > table_300test_${sample}_${sim}.RDPmidori_Arthropoda_order.txt
#           mv table_300test_${sample}_${sim}.RDPmidori_Arthropoda_order.txt table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt
#           echo "table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt"
#           grep "Stomatopoda" table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt
#      done
# done

# DANGEROUS CODE:  RUN ONLY ONCE AFTER GENERATING THE RDP ARTHROPODA-ONLY TABLES, BECAUSE IF RUN MORE THAN ONCE, WILL INSERT THE NEW TAXONOMIC RANK (e.g. Caliscelidae family 0.5)  MORE THAN ONCE
# for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
# do
#      for sim in `seq 97 97`
#      do
#           cd ${HOMEFOLDER}${ANALYSIS}/${OTUTABLEFOLDER}/OTU_tables
#           gsed -E 's/Frutillaria/Frutillaria\tgenus\t0.50\tFrutillaria/' table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt > table_300test_${sample}_${sim}.RDPmidori_Arthropoda_order.txt
#           mv table_300test_${sample}_${sim}.RDPmidori_Arthropoda_order.txt table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt
#           echo "table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt"
#           grep "Frutillaria" table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt
#      done
# done

# DANGEROUS CODE:  RUN ONLY ONCE AFTER GENERATING THE RDP ARTHROPODA-ONLY TABLES, BECAUSE IF RUN MORE THAN ONCE, WILL INSERT THE NEW TAXONOMIC RANK (e.g. Caliscelidae family 0.5)  MORE THAN ONCE
# for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
# do
#      for sim in `seq 97 97`
#      do
#           cd ${HOMEFOLDER}${ANALYSIS}/${OTUTABLEFOLDER}/OTU_tables
#           gsed -E 's/Cyrtanaspis\tgenus/Scraptiidae\tfamily\t0.50\tCyrtanaspis\tgenus/' table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt > table_300test_${sample}_${sim}.RDPmidori_Arthropoda_order.txt
#           mv table_300test_${sample}_${sim}.RDPmidori_Arthropoda_order.txt table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt
#           echo "table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt"
#           grep "Scraptiidae" table_300test_${sample}_${sim}.RDPmidori_Arthropoda.txt
#      done
# done
#### End script here.
exit

# Read into R and filter the OTUs via the phyloseq method and direct observation of the table. Also look for remove the Positive Control (PC) samples and any samples that are Illumina cross-talk.
# This is my final OTU table, with taxonomic assignments.

# Analyses:
# BLAST the OTUs against the MTB fasta to see drop outs and drop ins and to get insight into improvements in filtering.
# Procrustes analyses of all pairwise comparisons to test for tag bias

# https://github.com/hallamlab/mp_tutorial/wiki/Taxonomic-Analysis


# Use the list of retained OTUS to vsearch against the MTB ref dataset for analysis of dropouts and dropins
# Read into R and filter the OTUs via the phyloseq method and direct observation of the table. Also remove any samples that are Illumina cross-talk.
# This is my final OTU table, with taxonomic assignments.
# Analyses:

# BLAST the OTUs against the MTB fasta to see drop outs and drop ins and to get insight into improvements in filtering.
# Procrustes analyses of all pairwise comparisons to test for tag bias

# https://github.com/hallamlab/mp_tutorial/wiki/Taxonomic-Analysis

########################################
##### BLAST the final OTUs against the MTB reference sequences. Note that the MTB reference seqs are already clustered at 97% via sumaclust. There are 254 sumaclust 97% OTUs.

cd ${HOMEFOLDER}/data/MTB
# run once
# makeblastdb -in MTB_AllInputRefSeqs_20170726.fasta -dbtype nucl # to make the MTB ref dataset BLASTABLE

OTUTABLEFOLDER="OTUs_min2PCRs_min4copies_2017-11-26_time-1727"
echo $SUMASIM # should be 97

# take the top hit from blast output. visual inspection of the blast hits with similarity < 0.98 shows that the low-similarity hits are echo OTUs (echo OTUs are small OTUs that are similar to large OTUs that hit an MTB sequence at ~100% similarity)
parallel --jobs 4 "blastn -db ${HOMEFOLDER}data/MTB/MTB_AllInputRefSeqs_20170726.fasta -query ${HOMEFOLDER}analysis/${OTUTABLEFOLDER}/OTU_tables/table_300test_{1}_${SUMASIM}_Arthropoda.fas -num_threads 2 -evalue 1e-10 -max_target_seqs 1 -outfmt 6 -out ${HOMEFOLDER}analysis/${OTUTABLEFOLDER}/OTU_tables/table_300test_{1}_${SUMASIM}_Arthropoda.blastnMTB.txt" ::: A B C D E F

# single library code
# EXPERIMENT=F
# blastn -db ${HOMEFOLDER}/data/MTB/MTB_AllInputRefSeqs_20170726.fasta -query ${HOMEFOLDER}/analysis/${OTUTABLEFOLDER}/OTU_tables/table_300test_${EXPERIMENT}_${SUMASIM}_Arthropoda.fas -num_threads 6 -evalue 1e-10 -max_target_seqs 1 -outfmt 6 -out ${HOMEFOLDER}/analysis/${OTUTABLEFOLDER}/OTU_tables/table_300test_${EXPERIMENT}_${SUMASIM}_Arthropoda.blastnMTB.txt

# outfmt 6 column headings
# qseqid	sseqid	pident	length	mismatch	gapopen	qstart	qend	sstart	send	evalue	bitscore

# now run Dropout analyses in R

############################################################################################################
############ SINGLE POOL analysis #############
### OTU tables for single pool analysis
mkdir ${HOMEFOLDER}/analysis/singlepools/
cd ${HOMEFOLDER}data
# perl ${HOMEFOLDER}scripts/300Test_extractAllreadsfromSortV1.pl -id seqs/ -o sep_pools_v1.fas  # this version does not dereplicate the reads, so the file output is bigger
perl ${HOMEFOLDER}scripts/300Test_extractAllreadsfromSortV2.pl -id seqs/ -o sep_pools_v2.fas
# perl ${HOMEFOLDER}scripts/300Test_extractAllreadsfromSort.pl -h

vsearch --fastx_filter sep_pools_v2.fas --fastq_minlen ${MINLEN} --fastq_maxlen ${MAXLEN} --fastaout sep_pools_v2_${MINLEN}-${MAXLEN}.fas
usearch -fastx_info sep_pools_v2_${MINLEN}-${MAXLEN}.fas # confirm that the lengths are indeed limited
head sep_pools_v2.fas # to check format of header. Should look like this: The pool name is separated from the rest of the header (e.g. folderA1), and the index number is preceded by a :.  There should also be a count=NNNN at the end
# >folderF3 Tag9-Tag9:2662986 count=1
rm sep_pools_v2.fas

mv sep_pools_v2_${MINLEN}-${MAXLEN}.fas ${HOMEFOLDER}/analysis/singlepools/sep_pools_v2.fas
cd ${HOMEFOLDER}/analysis/singlepools/

# Make separate fasta files for each pool
sample_libs=($(cat ${HOMEFOLDER}data/seqs/samplelist.txt | cut -c 1 | uniq))  # cut out all but the first letter of each filename and keep unique values, samplelist.txt should already exist
echo "There are" "${#sample_libs[@]}" "samples that will be processed:"  "${sample_libs[@]}" # echo number and name of elements in the array
echo "There are ${POOLS} pools per library."

# make a separate file with the name of each pool
for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
              for pool in `seq 1 ${POOLS}`
              do
                            echo "folder${sample}${pool}" > folder${sample}${pool}.list
              done
done

# takes a minute
parallel -k 'seqtk subseq sep_pools_v2.fas folder{1}{2}.list > sep_pools_{1}{2}.fas' ::: "${sample_libs[@]}" ::: `seq 1 ${POOLS}`

# view header names to ensure that the correct pool is in each fasta file
for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
              for pool in `seq 1 ${POOLS}`
              do
                            echo "pool${sample}${pool}"
                            bioawk -c fastx '{ print ">"$name" "$comment }' sep_pools_${sample}${pool}.fas | head -n 1
                            bioawk -c fastx '{ print ">"$name" "$comment }' sep_pools_${sample}${pool}.fas | tail -n 1
              done
done

rm folder*.list
rm sep_pools_v2.fas

# to see more than 1 sequence per pool
# bioawk -c fastx '{ print ">"$name" "$comment"\n"$seq }' sep_pools_A1.fas.gz | head -n 3000 > A1.fas

# remove pool name (e.g. folderA1) from fasta headers
for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
              for pool in `seq 1 ${POOLS}`
              do
                            # gzip -d sep_pools_${sample}${pool}.fas.gz
                            gsed -E 's/>folder[A-F][1-3] />/' sep_pools_${sample}${pool}.fas > sep_pools_${sample}${pool}_folderremoved.fas
                            rm sep_pools_${sample}${pool}.fas
              done
done

# check header formats.  Should look like:  >Tag10-Tag10:2002610 count=190
for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
     for pool in `seq 1 ${POOLS}`
     do
          echo pool${sample}${pool}
          bioawk -c fastx '{ print ">"$name" "$comment }' sep_pools_${sample}${pool}_folderremoved.fas | head -n 1
          bioawk -c fastx '{ print ">"$name" "$comment }' sep_pools_${sample}${pool}_folderremoved.fas | tail -n 1
     done
done


# vsearch --uchime_denovo.  takes ~12-13 mins per file. I cannot run parallel on this, even where I run parallel on a file of commands
for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
         for pool in `seq 1 ${POOLS}`
         do
              gsed 's/ count/;size/' sep_pools_${sample}${pool}_folderremoved.fas > sep_pools_${sample}${pool}_forvsearch.fas;
              vsearch --sortbysize sep_pools_${sample}${pool}_forvsearch.fas --output sep_pools_${sample}${pool}_forvsearch_sorted.fas;
              vsearch --uchime_denovo sep_pools_${sample}${pool}_forvsearch_sorted.fas --nonchimeras sep_pools_${sample}${pool}_forvsearch_sorted_nochimeras.fas;
              gsed 's/;size/ count/' sep_pools_${sample}${pool}_forvsearch_sorted_nochimeras.fas > sep_pools_${sample}${pool}_forsumaclust.fas
         done
done

# # uchime2_denovo version - This takes > 6 hours per file to run on a large file like this
# rm uchime_commands.txt # ensure no uchime_commands.txt file is present
# # create a list of commands with the correct arguments
# for sample in ${sample_libs[@]}  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
# do
#          for pool in `seq 1 ${POOLS}`
#          do
#               gsed 's/ count/;size/' sep_pools_${sample}${pool}_folderremoved.fas > sep_pools_${sample}${pool}_foruchime.fas;
#               usearch9 -sortbysize sep_pools_${sample}${pool}_foruchime.fas -fastaout sep_pools_${sample}${pool}_foruchime_sorted.fas;
#               usearch9 -uchime2_denovo sep_pools_${sample}${pool}_foruchime_sorted.fas -nonchimeras sep_pools_${sample}${pool}_foruchime_sorted_nochimeras.fas;
#               gsed 's/;size/ count/' sep_pools_${sample}${pool}_foruchime_sorted_nochimeras.fas > sep_pools_${sample}${pool}_forsumaclust.fas
#          done
# done
# # run parallel --dryrun to see the commands that will be run, without actually running them.
# # cat uchime_commands.txt | wc -l # should be 72
# parallel -k --jobs 2 :::: uchime_commands.txt; rm uchime_commands.txt  # parallel :::: uchime_commands.txt means that the commands come from uchime_commands.txt.  I use --jobs 2 because hyperthreading might slow things down. After running, rm uchime_commands.txt.


# remove working files
for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
         for pool in `seq 1 ${POOLS}`
         do
              rm sep_pools_${sample}${pool}_forvsearch.fas
              rm sep_pools_${sample}${pool}_forvsearch_sorted.fas
              rm sep_pools_${sample}${pool}_forvsearch_sorted_nochimeras.fas
         done
done


# 96% sumaclust
echo ${SUMASIM} # confirm that there is a similarity value chosen
SUMASIM=96 # if there is no SUMASIM value
echo ${SUMASIM} # confirm that there is a similarity value chosen

# SUMACLUST
rm sumaclust_commands.txt # ensure no sumaclust_commands.txt file is present
# create a list of commands with the correct arguments
for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
         for pool in `seq 1 ${POOLS}`
         do
                  echo "${HOMEFOLDER}/analysis/singlepools/; \
                  sumaclust -t .${SUMASIM} -e sep_pools_${sample}${pool}_forsumaclust.fas > OTUs_${SUMASIM}_sumaclust_${sample}${pool}.fna; \
                  python ${DAME}tabulateSumaclust.py -i OTUs_${SUMASIM}_sumaclust_${sample}${pool}.fna -o table_300test_${SUMASIM}_${sample}${pool}.txt -blast; \
                  gzip OTUs_${SUMASIM}_sumaclust_${sample}${pool}.fna; \
                  gzip sep_pools_${sample}${pool}_forsumaclust.fas" \
                  >> sumaclust_commands.txt
         done
done

# change filename suffix of OTU seq files to *.fas
for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
         for pool in `seq 1 ${POOLS}`
         do
                  mv table_300test_${SUMASIM}_${sample}${pool}.txt.blast.txt table_300test_${SUMASIM}_${sample}${pool}.fas
         done
done

# 97% sumaclust
echo ${SUMASIM} # confirm that there is a similarity value chosen
SUMASIM=97 # if there is no SUMASIM value
echo ${SUMASIM} # confirm that there is a similarity value chosen

# SUMACLUST
rm sumaclust_commands.txt # ensure no sumaclust_commands.txt file is present
# create a list of commands with the correct arguments
for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
         for pool in `seq 1 ${POOLS}`
         do
                  echo "${HOMEFOLDER}/analysis/singlepools/; \
                  sumaclust -t .${SUMASIM} -e sep_pools_${sample}${pool}_forsumaclust.fas > OTUs_${SUMASIM}_sumaclust_${sample}${pool}.fna; \
                  python ${DAME}tabulateSumaclust.py -i OTUs_${SUMASIM}_sumaclust_${sample}${pool}.fna -o table_300test_${SUMASIM}_${sample}${pool}.txt -blast; \
                  gzip OTUs_${SUMASIM}_sumaclust_${sample}${pool}.fna; \
                  gzip sep_pools_${sample}${pool}_forsumaclust.fas" \
                  >> sumaclust_commands.txt
         done
done


# run parallel --dryrun to see the commands that will be run, without actually running them.
# cat sumaclust_commands.txt | wc -l # should be 18
parallel -k --jobs 4 :::: sumaclust_commands.txt  # parallel :::: sumaclust_commands.txt means that the commands come from sumaclust_commands.txt.  Can use --jobs 3 because my laptop's Intel i7 chip has 4 cores with 2 threads each.  Using 4 lets me use the remaining thread for other work.  After running, rm sumaclust_commands.txt to ensure that no command text remains.

# change filename suffix of OTU seq files to *.fas
for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
         for pool in `seq 1 ${POOLS}`
         do
                  mv table_300test_${SUMASIM}_${sample}${pool}.txt.blast.txt table_300test_${SUMASIM}_${sample}${pool}.fas
         done
done

rm sumaclust_commands.txt
rm sep_pools_*_folderremoved.fas



# LOOP VERSION
# for sample in ${sample_libs[@]}  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
# do
#               for pool in `seq 1 ${POOLS}`
#               do
#                             cd ${HOMEFOLDER}/analysis/singlepools/
#                             sumaclust -t .${SUMASIM} -e sep_pools_${sample}${pool}_folderremoved.fas > OTUs_${SUMASIM}_sumaclust_${sample}${pool}.fna
#                             python ${DAME}tabulateSumaclust.py -i OTUs_${SUMASIM}_sumaclust_${sample}${pool}.fna -o table_300test_${SUMASIM}_${sample}${pool}.txt -blast
#                             gzip OTUs_${SUMASIM}_sumaclust_${sample}${pool}.fna
#                             gzip sep_pools_${sample}${pool}_folderremoved.fas
#               done
# done

# BUGGY: PARALLEL COMPOSED COMMANDS VERSION:  parallel code with composed commands.  For some reason, the variable names (e.g. ${HOMEFOLDER}) come out as blanks in the commands. Variable names work in non-composed commands, but in composed commands, they don't work. So i have used a loop to create all the separate commands and use that as input in parallel
              # date; parallel -k 'cd ${HOMEFOLDER}/analysis/singlepools/; sumaclust -t .${SUMASIM} -e sep_pools_{1}{2}_folderremoved.fas > OTUs_${SUMASIM}_sumaclust_{1}{2}.fna; python ${DAME}tabulateSumaclust.py -i OTUs_${SUMASIM}_sumaclust_{1}{2}.fna -o table_300test_${SUMASIM}_{1}{2}.txt -blast' ::: ${sample_libs[@]} ::: `seq 1 ${POOLS}`; date

# # uchime2_denovo DON'T RUN THIS
# # runs very fast, but there's no size information to help uchime work
# for sample in ${sample_libs[@]}  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
# do
#          for pool in `seq 1 ${POOLS}`
#          do
#               gsed -E 's/(>OTU[0-9]+)/\1;size=1/' table_300test_96_${sample}${pool}.txt.blast.txt > table_300test_96_${sample}${pool}.size.blast.txt;
#               usearch9 -uchime2_denovo table_300test_96_${sample}${pool}.size.blast.txt -nonchimeras table_300test_96_${sample}${pool}.nochimeras.blast.txt;
#          done
# done


# After sumaclust, I upload the OTU fasta files to hpc.uea.ac.uk and RDP classify the sumaclust OTUs against the midori COI database, and I filter out all non-Arthropoda.  This has to run on the hpc.uea.ac.uk server because it needs a lot of RAM. The RDP Classifier was trained on the Midori COI database, and the trained data and reference sequences are on hpc:  ~/midori/

# to generate RDPClassify commands for the singlepool analysis
cd ${HOMEFOLDER}data/seqs
rm RDPClassify.txt
for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
              for pool in `seq 1 ${POOLS}`
              do
                            echo "java -Xmx64g -jar ~/scripts/RDPTools/classifier.jar classify -t ~/midori/RDP_out/rRNAClassifier.properties -o ~/Yahan/300Test/singlepools/table_300test_${sample}${pool}_${SUMASIM}.RDPmidori.txt ~/Yahan/300Test/singlepools/table_300test_${SUMASIM}_${sample}${pool}.fas" >> RDPClassify.txt
              done
done


# Upload to hpc.uea.ac.uk and make this bsub file
# make sure to upload to the correct folder

# requires ~10 mins per OTU fasta file
################ BSUB START ################
#!/bin/sh
#BSUB -q mellanox-ib     # on hpc.uea.ac.uk, mellanox-ib (168 hours = 7 days)
#BSUB -J RDPclmidAF
#BSUB -oo RDPclmidAF.out
#BSUB -eo RDPclmidAF.err
#BSUB -R "rusage[mem=64000]"  # max mem for mellanox-ib is 128GB = 128000
#BSUB -M 64000
#BSUB -B        # sends email to me when job starts
#BSUB -N        # sends email to me when job finishes
. /etc/profile
module load java/jdk1.8.0_51
java -Xmx64g -jar ~/scripts/RDPTools/classifier.jar classify -t ~/midori/RDP_out/rRNAClassifier.properties -o ~/Yahan/300Test/singlepools/table_300test_A1_97.RDPmidori.txt ~/Yahan/300Test/singlepools/table_300test_97_A1.fas
java -Xmx64g -jar ~/scripts/RDPTools/classifier.jar classify -t ~/midori/RDP_out/rRNAClassifier.properties -o ~/Yahan/300Test/singlepools/table_300test_A2_97.RDPmidori.txt ~/Yahan/300Test/singlepools/table_300test_97_A2.fas
java -Xmx64g -jar ~/scripts/RDPTools/classifier.jar classify -t ~/midori/RDP_out/rRNAClassifier.properties -o ~/Yahan/300Test/singlepools/table_300test_A3_97.RDPmidori.txt ~/Yahan/300Test/singlepools/table_300test_97_A3.fas
java -Xmx64g -jar ~/scripts/RDPTools/classifier.jar classify -t ~/midori/RDP_out/rRNAClassifier.properties -o ~/Yahan/300Test/singlepools/table_300test_B1_97.RDPmidori.txt ~/Yahan/300Test/singlepools/table_300test_97_B1.fas
java -Xmx64g -jar ~/scripts/RDPTools/classifier.jar classify -t ~/midori/RDP_out/rRNAClassifier.properties -o ~/Yahan/300Test/singlepools/table_300test_B2_97.RDPmidori.txt ~/Yahan/300Test/singlepools/table_300test_97_B2.fas
java -Xmx64g -jar ~/scripts/RDPTools/classifier.jar classify -t ~/midori/RDP_out/rRNAClassifier.properties -o ~/Yahan/300Test/singlepools/table_300test_B3_97.RDPmidori.txt ~/Yahan/300Test/singlepools/table_300test_97_B3.fas
java -Xmx64g -jar ~/scripts/RDPTools/classifier.jar classify -t ~/midori/RDP_out/rRNAClassifier.properties -o ~/Yahan/300Test/singlepools/table_300test_C1_97.RDPmidori.txt ~/Yahan/300Test/singlepools/table_300test_97_C1.fas
java -Xmx64g -jar ~/scripts/RDPTools/classifier.jar classify -t ~/midori/RDP_out/rRNAClassifier.properties -o ~/Yahan/300Test/singlepools/table_300test_C2_97.RDPmidori.txt ~/Yahan/300Test/singlepools/table_300test_97_C2.fas
java -Xmx64g -jar ~/scripts/RDPTools/classifier.jar classify -t ~/midori/RDP_out/rRNAClassifier.properties -o ~/Yahan/300Test/singlepools/table_300test_C3_97.RDPmidori.txt ~/Yahan/300Test/singlepools/table_300test_97_C3.fas
java -Xmx64g -jar ~/scripts/RDPTools/classifier.jar classify -t ~/midori/RDP_out/rRNAClassifier.properties -o ~/Yahan/300Test/singlepools/table_300test_D1_97.RDPmidori.txt ~/Yahan/300Test/singlepools/table_300test_97_D1.fas
java -Xmx64g -jar ~/scripts/RDPTools/classifier.jar classify -t ~/midori/RDP_out/rRNAClassifier.properties -o ~/Yahan/300Test/singlepools/table_300test_D2_97.RDPmidori.txt ~/Yahan/300Test/singlepools/table_300test_97_D2.fas
java -Xmx64g -jar ~/scripts/RDPTools/classifier.jar classify -t ~/midori/RDP_out/rRNAClassifier.properties -o ~/Yahan/300Test/singlepools/table_300test_D3_97.RDPmidori.txt ~/Yahan/300Test/singlepools/table_300test_97_D3.fas
java -Xmx64g -jar ~/scripts/RDPTools/classifier.jar classify -t ~/midori/RDP_out/rRNAClassifier.properties -o ~/Yahan/300Test/singlepools/table_300test_E1_97.RDPmidori.txt ~/Yahan/300Test/singlepools/table_300test_97_E1.fas
java -Xmx64g -jar ~/scripts/RDPTools/classifier.jar classify -t ~/midori/RDP_out/rRNAClassifier.properties -o ~/Yahan/300Test/singlepools/table_300test_E2_97.RDPmidori.txt ~/Yahan/300Test/singlepools/table_300test_97_E2.fas
java -Xmx64g -jar ~/scripts/RDPTools/classifier.jar classify -t ~/midori/RDP_out/rRNAClassifier.properties -o ~/Yahan/300Test/singlepools/table_300test_E3_97.RDPmidori.txt ~/Yahan/300Test/singlepools/table_300test_97_E3.fas
java -Xmx64g -jar ~/scripts/RDPTools/classifier.jar classify -t ~/midori/RDP_out/rRNAClassifier.properties -o ~/Yahan/300Test/singlepools/table_300test_F1_97.RDPmidori.txt ~/Yahan/300Test/singlepools/table_300test_97_F1.fas
java -Xmx64g -jar ~/scripts/RDPTools/classifier.jar classify -t ~/midori/RDP_out/rRNAClassifier.properties -o ~/Yahan/300Test/singlepools/table_300test_F2_97.RDPmidori.txt ~/Yahan/300Test/singlepools/table_300test_97_F2.fas
java -Xmx64g -jar ~/scripts/RDPTools/classifier.jar classify -t ~/midori/RDP_out/rRNAClassifier.properties -o ~/Yahan/300Test/singlepools/table_300test_F3_97.RDPmidori.txt ~/Yahan/300Test/singlepools/table_300test_97_F3.fas
# add all commands from RDPClassify.txt and run on hpc
################ BSUB END ################

cd ${HOMEFOLDER}data/seqs
rm RDPClassify.txt

# 9.
# Download files to mac: e.g. table_300test_A1_97.RDPmidori.txt
# Filter out non-Arthropoda from RDP assignment table.  Keep only Arthropoda with prob >= ARTHMINPROB (set to 0.80).
# Filter out non-Arthropoda OTUs from OTU representative sequences fasta file


OTUTABLEFOLDER="singlepools"
cd ${HOMEFOLDER}${ANALYSIS}/${OTUTABLEFOLDER}
for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
     for sim in `seq 97 97`  # change to `seq 96 97` if i want to do both 96 and 97
     do
          for pool in `seq 1 ${POOLS}`
          do
          awk -v arthmin=${ARTHMINPROB} '$8 ~ /Arthropoda/ && $10 >= arthmin { print }' table_300test_${sample}${pool}_${sim}.RDPmidori.txt > table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt
          gsed -E 's/\t\t/\t/' table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt > table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda_nodbltab.txt  # remove double tab after first column
          mv table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda_nodbltab.txt table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt
          seqtk subseq table_300test_${sim}_${sample}${pool}.fas <(cut -f 1 table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt) > table_300test_${sim}_${sample}${pool}_Arthropoda.fas # filter the OTU representative sequences file (e.g. table_300test_97_A1.fas) to include only the Arthropoda. Output fasta file:  table_300test_97_A1_Arthropoda.fas
               # <(cut -f 1 table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt) produces the list of sequences to keep
          done
     done
done

# Check that the right number of OTUs has been removed from each file
for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
     for sim in `seq 97 97`  # so far, have only done 96% similarity
     do
          for pool in `seq 1 ${POOLS}`
          do
          echo "OTU tables"
          wc -l table_300test_${sample}${pool}_${sim}.RDPmidori.txt
          wc -l table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt
          echo "fasta files"
          grep ">" table_300test_${sim}_${sample}${pool}.fas | wc -l
          grep ">" table_300test_${sim}_${sample}${pool}_Arthropoda.fas | wc -l
          done
     done
done

#### Prepare OTU tables for tag-bias analysis in R:  300Test_singlepools_OTU_table.R.  This code is to replace the headings in the singlepool OTU tables from TagN-TagN to the actual sample names as listed in PSinfo. However, The Midori database has a problem because NCBI has a problem.  Some of the taxonomic assignments are not complete. For example, the Bruchomorpha records do not have the family name:  Caliscelidae. This creates an OTU table with fewer columns (e.g. missing:  Caliscelidae   family    0.n). Thus, i use this code to add this missing bits.

# This is a bit complicated. First, I try to import the rdpmidori files into R, but some of them will fail because the number of fields is incorrect (i.e. not 25 fields).  So i look at the row number of the file that causes R to crash and write/run code to fix that (see below)

# DANGEROUS CODE:  RUN ONLY ONCE AFTER GENERATING THE RDP ARTHROPODA-ONLY TABLES, BECAUSE IF RUN MORE THAN ONCE, WILL INSERT THE NEW TAXONOMIC RANK (e.g. Caliscelidae family 0.5) MORE THAN ONCE
# OTUTABLEFOLDER="singlepools"
# for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
# do
#      for sim in `seq 97 97`  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
#      do
#           for pool in `seq 1 ${POOLS}`
#           do
#           cd ${HOMEFOLDER}${ANALYSIS}/${OTUTABLEFOLDER}
#           gsed -E 's/Bruchomorpha\tgenus/Caliscelidae\tfamily\t0.50\t\Bruchomorpha\tgenus/' table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt > table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda_family.txt
#           mv table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda_family.txt table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt
#           echo "table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt"
#           grep "Bruchomorpha" table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt
#           done
#      done
# done

# DANGEROUS CODE:  RUN ONLY ONCE AFTER GENERATING THE RDP ARTHROPODA-ONLY TABLES, BECAUSE IF RUN MORE THAN ONCE, WILL INSERT THE NEW TAXONOMIC RANK (e.g. Caliscelidae family 0.5)  MORE THAN ONCE
# OTUTABLEFOLDER="singlepools"
# for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
# do
#      for sim in `seq 97 97`
#      do
#           for pool in `seq 1 ${POOLS}`
#           do
#           cd ${HOMEFOLDER}${ANALYSIS}/${OTUTABLEFOLDER}
#           gsed -E 's/Hypsidia\tgenus/Drepanidae\tfamily\t0.50\tHypsidia\tgenus/' table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt > table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda_order.txt
#           mv table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda_order.txt table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt
#           echo "table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt"
#           grep "Drepanidae" table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt
#           done
#      done
# done

# DANGEROUS CODE:  RUN ONLY ONCE AFTER GENERATING THE RDP ARTHROPODA-ONLY TABLES, BECAUSE IF RUN MORE THAN ONCE, WILL INSERT THE NEW TAXONOMIC RANK (e.g. Caliscelidae family 0.5)  MORE THAN ONCE
# OTUTABLEFOLDER="singlepools"
# for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
# do
#      for sim in `seq 97 97`
#      do
#           for pool in `seq 1 ${POOLS}`
#           do
#           cd ${HOMEFOLDER}${ANALYSIS}/${OTUTABLEFOLDER}
#           gsed -E 's/Choroterpes\tgenus/Leptophlebiidae\tfamily\t0.50\tChoroterpes\tgenus/' table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt > table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda_order.txt
#           mv table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda_order.txt table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt
#           echo "table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt"
#           grep "Leptophlebiidae" table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt
#           done
#      done
# done

# DANGEROUS CODE:  RUN ONLY ONCE AFTER GENERATING THE RDP ARTHROPODA-ONLY TABLES, BECAUSE IF RUN MORE THAN ONCE, WILL INSERT THE NEW TAXONOMIC RANK (e.g. Caliscelidae family 0.5)  MORE THAN ONCE
# OTUTABLEFOLDER="singlepools"
# for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
# do
#      for sim in `seq 97 97`
#      do
#           for pool in `seq 1 ${POOLS}`
#           do
#           cd ${HOMEFOLDER}${ANALYSIS}/${OTUTABLEFOLDER}
#           gsed -E 's/Plectrotarsus\tgenus/Plectrotarsidae\tfamily\t0.50\tPlectrotarsus\tgenus/' table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt > table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda_order.txt
#           mv table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda_order.txt table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt
#           echo "table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt"
#           grep "Plectrotarsidae" table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt
#           done
#      done
# done

# DANGEROUS CODE:  RUN ONLY ONCE AFTER GENERATING THE RDP ARTHROPODA-ONLY TABLES, BECAUSE IF RUN MORE THAN ONCE, WILL INSERT THE NEW TAXONOMIC RANK (e.g. Caliscelidae family 0.5)  MORE THAN ONCE
# OTUTABLEFOLDER="singlepools"
# for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
# do
#      for sim in `seq 97 97`
#      do
#           for pool in `seq 1 ${POOLS}`
#           do
#           cd ${HOMEFOLDER}${ANALYSIS}/${OTUTABLEFOLDER}
#           gsed -E 's/Lepidotrichidae\tfamily/Zygentoma\torder\t0.50\tLepidotrichidae\tfamily/' table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt > table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda_order.txt
#           mv table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda_order.txt table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt
#           echo "table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt"
#           grep "Lepidotrichidae" table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt
#           done
#      done
# done

# DANGEROUS CODE:  RUN ONLY ONCE AFTER GENERATING THE RDP ARTHROPODA-ONLY TABLES, BECAUSE IF RUN MORE THAN ONCE, WILL INSERT THE NEW TAXONOMIC RANK (e.g. Caliscelidae family 0.5)  MORE THAN ONCE
# OTUTABLEFOLDER="singlepools"
# for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
# do
#      for sim in `seq 97 97`
#      do
#           for pool in `seq 1 ${POOLS}`
#           do
#           cd ${HOMEFOLDER}${ANALYSIS}/${OTUTABLEFOLDER}
#           gsed -E 's/Entomobryidae\tfamily/Entomobryomorpha\torder\t0.50\tEntomobryidae\tfamily/' table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt > table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda_order.txt
#           mv table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda_order.txt table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt
#           echo "table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt"
#           grep "Entomobryomorpha" table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt
#           done
#      done
# done

# DANGEROUS CODE:  RUN ONLY ONCE AFTER GENERATING THE RDP ARTHROPODA-ONLY TABLES, BECAUSE IF RUN MORE THAN ONCE, WILL INSERT THE NEW TAXONOMIC RANK (e.g. Caliscelidae family 0.5)  MORE THAN ONCE
# OTUTABLEFOLDER="singlepools"
# for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
# do
#      for sim in `seq 97 97`
#      do
#           for pool in `seq 1 ${POOLS}`
#           do
#           cd ${HOMEFOLDER}${ANALYSIS}/${OTUTABLEFOLDER}
#           gsed -E 's/Habrosyne\tgenus/Drepanidae\tfamily\t0.50\tHabrosyne\tgenus/' table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt > table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda_order.txt
#           mv table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda_order.txt table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt
#           echo "table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt"
#           grep "Drepanidae" table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt
#           done
#      done
# done

# DANGEROUS CODE:  RUN ONLY ONCE AFTER GENERATING THE RDP ARTHROPODA-ONLY TABLES, BECAUSE IF RUN MORE THAN ONCE, WILL INSERT THE NEW TAXONOMIC RANK (e.g. Caliscelidae family 0.5)  MORE THAN ONCE
# OTUTABLEFOLDER="singlepools"
# for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
# do
#      for sim in `seq 97 97`
#      do
#           for pool in `seq 1 ${POOLS}`
#           do
#           cd ${HOMEFOLDER}${ANALYSIS}/${OTUTABLEFOLDER}
#           gsed -E 's/Coronididae\tfamily/Stomatopoda\torder\t0.50\tCoronididae\tfamily/' table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt > table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda_order.txt
#           mv table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda_order.txt table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt
#           echo "table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt"
#           grep "Stomatopoda" table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt
#           done
#      done
# done

# DANGEROUS CODE:  RUN ONLY ONCE AFTER GENERATING THE RDP ARTHROPODA-ONLY TABLES, BECAUSE IF RUN MORE THAN ONCE, WILL INSERT THE NEW TAXONOMIC RANK (e.g. Caliscelidae family 0.5)  MORE THAN ONCE
# OTUTABLEFOLDER="singlepools"
# for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
# do
#      for sim in `seq 97 97`
#      do
#           for pool in `seq 1 ${POOLS}`
#           do
#           cd ${HOMEFOLDER}${ANALYSIS}/${OTUTABLEFOLDER}
#           gsed -E 's/Frutillaria/Frutillaria\tgenus\t0.50\tFrutillaria/' table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt > table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda_order.txt
#           mv table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda_order.txt table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt
#           echo "table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt"
#           grep "Frutillaria" table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt
#           done
#      done
# done

# DANGEROUS CODE:  RUN ONLY ONCE AFTER GENERATING THE RDP ARTHROPODA-ONLY TABLES, BECAUSE IF RUN MORE THAN ONCE, WILL INSERT THE NEW TAXONOMIC RANK (e.g. Caliscelidae family 0.5)  MORE THAN ONCE
# OTUTABLEFOLDER="singlepools"
# for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
# do
#      for sim in `seq 97 97`
#      do
#           for pool in `seq 1 ${POOLS}`
#           do
#           cd ${HOMEFOLDER}${ANALYSIS}/${OTUTABLEFOLDER}
#           gsed -E 's/Cyrtanaspis\tgenus/Scraptiidae\tfamily\t0.50\tCyrtanaspis\tgenus/' table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt > table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda_order.txt
#           mv table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda_order.txt table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt
#           echo "table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt"
#           grep "Scraptiidae" table_300test_${sample}${pool}_${sim}.RDPmidori_Arthropoda.txt
#           done
#      done
# done

#### Prepare OTU tables for tag-bias analysis in R:  300Test_singlepools_OTU_table.R.  This code is to replace the headings in the singlepool OTU tables from TagN-TagN to the actual sample names as listed in PSinfo
# e.g. I change the Tag2-Tag2 column heading changed to Hhmlbody. One complication is that different libraries use different tag pairs to indicate the same sample, but all this information is recorded in the PSinfo files.

# Step 1. For safety, run this code in a separate folder, copying in the files and then copying out the outputs:  tag_to_samples/.
# mkdir ${HOMEFOLDER}analysis/tags_to_samples
cd ${HOMEFOLDER}analysis/tags_to_samples/

# Step 2. copy the PSinfo_300test_COI files from ${HOMEFOLDER}data to ${HOMEFOLDER}analysis/tag_to_samples/, sort PSinfo files by pool number, and change the Tag format from "TagN TagN" to "TagN-TagN". Store in a new file.  The sorting isn't really necessary.
cp ${HOMEFOLDER}data/PSinfo_300test_COI*.txt ${HOMEFOLDER}analysis/tags_to_samples

sort -k4 PSinfo_300test_COIA.txt > PSinfo_300test_COIA_sorted.txt # sort by pool
gsed -E 's/(Tag[^\t]+)\tTag/\1-Tag/' PSinfo_300test_COIA_sorted.txt > PSinfo_300test_COIA_sorted_gsed.txt # put dash between Tagnames

sort -k4 PSinfo_300test_COIB.txt > PSinfo_300test_COIB_sorted.txt # sort by pool
gsed -E 's/(Tag[^\t]+)\tTag/\1-Tag/' PSinfo_300test_COIB_sorted.txt > PSinfo_300test_COIB_sorted_gsed.txt # put dash between Tagnames

sort -k4 PSinfo_300test_COIC.txt > PSinfo_300test_COIC_sorted.txt # sort by pool
gsed -E 's/(Tag[^\t]+)\tTag/\1-Tag/' PSinfo_300test_COIC_sorted.txt > PSinfo_300test_COIC_sorted_gsed.txt # put dash between Tagnames

sort -k4 PSinfo_300test_COID.txt > PSinfo_300test_COID_sorted.txt # sort by pool
gsed -E 's/(Tag[^\t]+)\tTag/\1-Tag/' PSinfo_300test_COID_sorted.txt > PSinfo_300test_COID_sorted_gsed.txt # put dash between Tagnames

sort -k4 PSinfo_300test_COIE.txt > PSinfo_300test_COIE_sorted.txt # sort by pool
gsed -E 's/(Tag[^\t]+)\tTag/\1-Tag/' PSinfo_300test_COIE_sorted.txt > PSinfo_300test_COIE_sorted_gsed.txt # put dash between Tagnames

sort -k4 PSinfo_300test_COIF.txt > PSinfo_300test_COIF_sorted.txt # sort by pool
gsed -E 's/(Tag[^\t]+)\tTag/\1-Tag/' PSinfo_300test_COIF_sorted.txt > PSinfo_300test_COIF_sorted_gsed.txt # put dash between Tagnames

rm PSinfo_300test_COI{A,B,C,D,E,F}_sorted.txt # rm working files
rm PSinfo_300test_COI{A,B,C,D,E,F}.txt # rm copies of the PSinfo_300test_COI files

# Step 3. Split the PSinfo files into three parts, one for each pool, named like this:  PSinfo_300test_COIA_sorted_gsed_1.txt
# There are now 18 PSinfo files
sample_libs=($(cat ${HOMEFOLDER}data/seqs/samplelist.txt | cut -c 1 | uniq))  # cut out all but the first letter of each filename and keep unique values, samplelist.txt should already exist
echo "There are" "${#sample_libs[@]}" "samples that will be processed:"  "${sample_libs[@]}" # echo number and name of elements in the array
echo "There are ${POOLS} pools per library."

for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
     for pool in `seq 1 ${POOLS}`
     do
          awk -v poolnum=${pool} '$3 ~ poolnum { print $0 }' PSinfo_300test_COI${sample}_sorted_gsed.txt > PSinfo_300test_COI${sample}_sorted_gsed_${pool}.txt # match column 3 to pool number and print to new subfile
          lines=$(cat PSinfo_300test_COI${sample}_sorted_gsed_${pool}.txt | wc -l)  # should always be 13 in 300Test
          assert_eq ${lines} 13  # this uses the assert.sh code that i sourced at the top of the script
          if [ "$?" == 0 ]; then
               echo "PSinfo_300test_COI${sample}_sorted_gsed_${pool}.txt has 13 lines:  good"
          else
               echo "PSinfo_300test_COI${sample}_sorted_gsed_${pool}.txt does not have 13 lines:  bad"
          fi
     done
done

# Step 4. change headings on OTU table from TagN-TagN format to actual sample names as listed in the PSinfo files

cp ${HOMEFOLDER}analysis/singlepools/table_300test_{96,97}_{A,B,C,D,E,F}{1,2,3}.txt ${HOMEFOLDER}analysis/tags_to_samples

# This takes ~ 1 minute (surprisingly), because there is a lot of file saving, one for each substitution: 13*18*2=468, plus 468 file backups saved
for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
     for sim in `seq 97 97`
     do
          for pool in `seq 1 ${POOLS}`
          do
               for line in `seq 1 13`
               do
                    tagpair=$(cat PSinfo_300test_COI${sample}_sorted_gsed_${pool}.txt | head -${line} | tail -1 | cut -f 2) # "head -n | tail-1 | cut -f m" to extract mth field in ith line, but note that DANGER is that if use a value of n beyond last line, it will return value from last line
                    samplename=$(cat PSinfo_300test_COI${sample}_sorted_gsed_${pool}.txt | head -${line} | tail -1 | cut -f 1)
                    sed -i ".tmp" "s/${tagpair}/${samplename}/" table_300test_${sim}_${sample}${pool}.txt # sed -i ".tmp" is for changing a file in place and also making a backup file with ".tmp" as the suffix.
               done
          done
     done
done

rm table_300test_*.tmp # rm the sed backup files

# Step 5. change filename of OTU tables to reflect that samplenames were added. I do this as a separate step because I ran sed -i 13 times on the same file, and I need to finish all that before changing the name of the output file.
for sample in "${sample_libs[@]}"  # ${sample_libs[@]} is the full bash array: A,B,C,D,E,F.  So loop over all samples
do
     for sim in `seq 97 97`
     do
          for pool in `seq 1 ${POOLS}`
          do
               mv table_300test_${sim}_${sample}${pool}.txt table_300test_${sim}_${sample}${pool}_samplenames.txt
          done
     done
done

# copy output files back to singlepools/ folder
cp ${HOMEFOLDER}analysis/tags_to_samples/table_300test_{96,97}_{A,B,C,D,E,F}{1,2,3}_samplenames.txt ${HOMEFOLDER}analysis/singlepools/

# empty tags_to_samples/ folder
rm PSinfo_300test_COI{A,B,C,D,E,F}_sorted_gsed_{1,2,3}.txt # to prevent problems for next time
rm PSinfo_300test_COI{A,B,C,D,E,F}_sorted_gsed.txt # to prevent problems for next time
rm table_300test_{96,97}_{A,B,C,D,E,F}{1,2,3}_samplenames.txt

# You're now ready to go with the R scripts to filter the OTU tables by taxonomy, presence in Positive control samples, contamination as revealed by extraction blanks, and by min-OTU size, as analysed by phyloseq.
