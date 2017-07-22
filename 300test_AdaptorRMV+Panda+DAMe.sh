#!/bin/bash
set -e
set -u
set -o pipefail
#######################################################################################
#######################################################################################
# a shell script to loop through a set of amplicon files, denoise, and apply the DAMe protocol
#######################################################################################
#######################################################################################

# Usage: bash 300test_AdaptorRMV+Panda+DAMe.sh
# e.g. 300test_AdaptorRMV+Panda+DAMe.sh

PIPESTART=$(date)

# run script from ~/PATH_TO/300Test/scripts
# fastq files in 300Test/data/seqs/
# DAMe files in 300Test/data/
# analysis outputs in 300Test/analysis/
# this script in 300Test/scripts/

# set variables
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
echo "There are" ${#sample_names[@]} "samples that will be processed." # echo number of elements in the array

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
              sample_prefix="$( basename $sample "_L001_R1_001.fastq.gz")"
              echo ${sample_prefix} >> sickle.out
              echo ${sample_prefix}
              sickle pe -f Adaptermv${sample_prefix}.pair1.truncated -r Adaptermv${sample_prefix}.pair2.truncated -o sickle_${sample_prefix}_R1.fq -p sickle_${sample_prefix}_R2.fq -t sanger -s sickle_Single${sample_prefix}.fq >> sickle.out
done

##-s sickle_Single${sample_prefix}.fq:  keep single read （not paired but good sequences), will not be used in DAMe pipeline as we are using twin tags but can be used in genome research.

#### remove the AdapterRemoval files
# rm Adaptermv*.*


#### SPAdes for error correction, via BayesHammer
#### to do:  test bfc
#### choosing not to use spades.py --meta because amplicon data
# download from http://cab.spbu.ru/software/spades/
    # sudo cp -r SPAdes-3.10.1 /usr/local/bin
# or
    # brew install homebrew/science/spades

# loop over all samples and apply BayesHammer error correction
for sample in ${sample_names[@]}  # ${sample_names[@]} is the full bash array.  So loop over all samples
do
              sample_prefix="$( basename $sample "_L001_R1_001.fastq.gz")"
              echo ${sample_prefix}
              spades.py --only-error-correction -1 sickle_${sample_prefix}_R1.fq -2 sickle_${sample_prefix}_R2.fq -s sickle_Single${sample_prefix}.fq -o SPAdes_hammer${sample_prefix}
done

# START HERE:  check if spades.py has finished for all samples. If so, delete the spades commands below and start with pandaseq.
# Also,




# A1 to C3 can be pasted at the command line in a single batch
date;spades.py --only-error-correction -1 sickle_A1R1.fq -2 sickle_A1R2.fq -s sickle_SingleA1.fq -o SPAdes_hammerA1 ;date
date;spades.py --only-error-correction -1 sickle_A2R1.fq -2 sickle_A2R2.fq -s sickle_SingleA2.fq -o SPAdes_hammerA2 ;date
date;spades.py --only-error-correction -1 sickle_A3R1.fq -2 sickle_A3R2.fq -s sickle_SingleA3.fq -o SPAdes_hammerA3 ;date
date;spades.py --only-error-correction -1 sickle_B1R1.fq -2 sickle_B1R2.fq -s sickle_SingleB1.fq -o SPAdes_hammerB1 ;date
date;spades.py --only-error-correction -1 sickle_B2R1.fq -2 sickle_B2R2.fq -s sickle_SingleB2.fq -o SPAdes_hammerB2 ;date
date;spades.py --only-error-correction -1 sickle_B3R1.fq -2 sickle_B3R2.fq -s sickle_SingleB3.fq -o SPAdes_hammerB3 ;date
date;spades.py --only-error-correction -1 sickle_C1R1.fq -2 sickle_C1R2.fq -s sickle_SingleC1.fq -o SPAdes_hammerC1 ;date
date;spades.py --only-error-correction -1 sickle_C2R1.fq -2 sickle_C2R2.fq -s sickle_SingleC2.fq -o SPAdes_hammerC2 ;date
date;spades.py --only-error-correction -1 sickle_C3R1.fq -2 sickle_C3R2.fq -s sickle_SingleC3.fq -o SPAdes_hammerC3 ;date
# D1 to F3 can be pasted at the command line in a single batch
date;spades.py --only-error-correction -1 sickle_D1R1.fq -2 sickle_D1R2.fq -s sickle_SingleD1.fq -o SPAdes_hammerD1 ;date
date;spades.py --only-error-correction -1 sickle_D2R1.fq -2 sickle_D2R2.fq -s sickle_SingleD2.fq -o SPAdes_hammerD2 ;date
date;spades.py --only-error-correction -1 sickle_D3R1.fq -2 sickle_D3R2.fq -s sickle_SingleD3.fq -o SPAdes_hammerD3 ;date
date;spades.py --only-error-correction -1 sickle_E1R1.fq -2 sickle_E1R2.fq -s sickle_SingleE1.fq -o SPAdes_hammerE1 ;date
date;spades.py --only-error-correction -1 sickle_E2R1.fq -2 sickle_E2R2.fq -s sickle_SingleE2.fq -o SPAdes_hammerE2 ;date
date;spades.py --only-error-correction -1 sickle_E3R1.fq -2 sickle_E3R2.fq -s sickle_SingleE3.fq -o SPAdes_hammerE3 ;date
date;spades.py --only-error-correction -1 sickle_F1R1.fq -2 sickle_F1R2.fq -s sickle_SingleF1.fq -o SPAdes_hammerF1 ;date
date;spades.py --only-error-correction -1 sickle_F2R1.fq -2 sickle_F2R2.fq -s sickle_SingleF2.fq -o SPAdes_hammerF2 ;date
date;spades.py --only-error-correction -1 sickle_F3R1.fq -2 sickle_F3R2.fq -s sickle_SingleF3.fq -o SPAdes_hammerF3 ;date

#### remove the sickle files
# rm sickle_*.fq

#### PandaSeq  use for read overlapping
# download from https://github.com/neufeld/pandaseq/releases
    # PANDAseq-2.11.pkg
# or
    # brew install homebrew/science/pandaseq

for sample in ${sample_names[@]}  # ${sample_names[@]} is the full bash array.  So loop over all samples
do
    sample_prefix="$( basename $sample "_L001_R1_001.fastq.gz")"
    echo ${sample_prefix}
    pandaseq -f SPAdes_hammer_${sample_prefix}/corrected/sickle_${sample_prefix}_R1.00.0_0.cor.fastq.gz -r SPAdes_hammer_${sample_prefix}/corrected/sickle_${sample_prefix}_R2.00.0_0.cor.fastq.gz -A simple_bayesian -F -N -T 7 -g pandaseq_log_${sample_prefix}.txt -w sickle_cor_panda_${sample_prefix}.fastq
done


pandaseq -f SPAdes_hammerA1/corrected/sickle_A1R1.00.0_0.cor.fastq.gz -r SPAdes_hammerA1/corrected/sickle_A1R2.00.0_0.cor.fastq.gz -A simple_bayesian -F -N -T 7 -g pandaseq_logA1.txt -w sickle_cor_panda_A1.fastq
pandaseq -f SPAdes_hammerA2/corrected/sickle_A2R1.00.0_0.cor.fastq.gz -r SPAdes_hammerA2/corrected/sickle_A2R2.00.0_0.cor.fastq.gz -A simple_bayesian -F -N -T 7 -g pandaseq_logA2.txt -w sickle_cor_panda_A2.fastq
pandaseq -f SPAdes_hammerA3/corrected/sickle_A3R1.00.0_0.cor.fastq.gz -r SPAdes_hammerA3/corrected/sickle_A3R2.00.0_0.cor.fastq.gz -A simple_bayesian -F -N -T 7 -g pandaseq_logA3.txt -w sickle_cor_panda_A3.fastq
pandaseq -f SPAdes_hammerB1/corrected/sickle_B1R1.00.0_0.cor.fastq.gz -r SPAdes_hammerB1/corrected/sickle_B1R2.00.0_0.cor.fastq.gz -A simple_bayesian -F -N -T 7 -g pandaseq_logB1.txt -w sickle_cor_panda_B1.fastq
pandaseq -f SPAdes_hammerB2/corrected/sickle_B2R1.00.0_0.cor.fastq.gz -r SPAdes_hammerB2/corrected/sickle_B2R2.00.0_0.cor.fastq.gz -A simple_bayesian -F -N -T 7 -g pandaseq_logB2.txt -w sickle_cor_panda_B2.fastq
pandaseq -f SPAdes_hammerB3/corrected/sickle_B3R1.00.0_0.cor.fastq.gz -r SPAdes_hammerB3/corrected/sickle_B3R2.00.0_0.cor.fastq.gz -A simple_bayesian -F -N -T 7 -g pandaseq_logB3.txt -w sickle_cor_panda_B3.fastq
pandaseq -f SPAdes_hammerC1/corrected/sickle_C1R1.00.0_0.cor.fastq.gz -r SPAdes_hammerC1/corrected/sickle_C1R2.00.0_0.cor.fastq.gz -A simple_bayesian -F -N -T 7 -g pandaseq_logC1.txt -w sickle_cor_panda_C1.fastq
pandaseq -f SPAdes_hammerC2/corrected/sickle_C2R1.00.0_0.cor.fastq.gz -r SPAdes_hammerC2/corrected/sickle_C2R2.00.0_0.cor.fastq.gz -A simple_bayesian -F -N -T 7 -g pandaseq_logC2.txt -w sickle_cor_panda_C2.fastq
pandaseq -f SPAdes_hammerC3/corrected/sickle_C3R1.00.0_0.cor.fastq.gz -r SPAdes_hammerC3/corrected/sickle_C3R2.00.0_0.cor.fastq.gz -A simple_bayesian -F -N -T 7 -g pandaseq_logC3.txt -w sickle_cor_panda_C3.fastq
pandaseq -f SPAdes_hammerD1/corrected/sickle_D1R1.00.0_0.cor.fastq.gz -r SPAdes_hammerD1/corrected/sickle_D1R2.00.0_0.cor.fastq.gz -A simple_bayesian -F -N -T 7 -g pandaseq_logD1.txt -w sickle_cor_panda_D1.fastq
pandaseq -f SPAdes_hammerD2/corrected/sickle_D2R1.00.0_0.cor.fastq.gz -r SPAdes_hammerD2/corrected/sickle_D2R2.00.0_0.cor.fastq.gz -A simple_bayesian -F -N -T 7 -g pandaseq_logD2.txt -w sickle_cor_panda_D2.fastq
pandaseq -f SPAdes_hammerD3/corrected/sickle_D3R1.00.0_0.cor.fastq.gz -r SPAdes_hammerD3/corrected/sickle_D3R2.00.0_0.cor.fastq.gz -A simple_bayesian -F -N -T 7 -g pandaseq_logD3.txt -w sickle_cor_panda_D3.fastq
pandaseq -f SPAdes_hammerE1/corrected/sickle_E1R1.00.0_0.cor.fastq.gz -r SPAdes_hammerE1/corrected/sickle_E1R2.00.0_0.cor.fastq.gz -A simple_bayesian -F -N -T 7 -g pandaseq_logE1.txt -w sickle_cor_panda_E1.fastq
pandaseq -f SPAdes_hammerE2/corrected/sickle_E2R1.00.0_0.cor.fastq.gz -r SPAdes_hammerE2/corrected/sickle_E2R2.00.0_0.cor.fastq.gz -A simple_bayesian -F -N -T 7 -g pandaseq_logE2.txt -w sickle_cor_panda_E2.fastq
pandaseq -f SPAdes_hammerE3/corrected/sickle_E3R1.00.0_0.cor.fastq.gz -r SPAdes_hammerE3/corrected/sickle_E3R2.00.0_0.cor.fastq.gz -A simple_bayesian -F -N -T 7 -g pandaseq_logE3.txt -w sickle_cor_panda_E3.fastq
pandaseq -f SPAdes_hammerF1/corrected/sickle_F1R1.00.0_0.cor.fastq.gz -r SPAdes_hammerF1/corrected/sickle_F1R2.00.0_0.cor.fastq.gz -A simple_bayesian -F -N -T 7 -g pandaseq_logF1.txt -w sickle_cor_panda_F1.fastq
pandaseq -f SPAdes_hammerF2/corrected/sickle_F2R1.00.0_0.cor.fastq.gz -r SPAdes_hammerF2/corrected/sickle_F2R2.00.0_0.cor.fastq.gz -A simple_bayesian -F -N -T 7 -g pandaseq_logF2.txt -w sickle_cor_panda_F2.fastq
pandaseq -f SPAdes_hammerF3/corrected/sickle_F3R1.00.0_0.cor.fastq.gz -r SPAdes_hammerF3/corrected/sickle_F3R2.00.0_0.cor.fastq.gz -A simple_bayesian -F -N -T 7 -g pandaseq_logF3.txt -w sickle_cor_panda_F3.fastq

# remove SPAdes output
# rm -rf SPAdes_hammer*/
# compress pandseq_log
# gzip pandaseq_log*.txt

# gzip sickle_cor_panda_{A,B,C,D,E,F}{1,2,3}.fastq

gzip sickle_cor_panda_A1.fastq
gzip sickle_cor_panda_A2.fastq
gzip sickle_cor_panda_A3.fastq
gzip sickle_cor_panda_B1.fastq
gzip sickle_cor_panda_B2.fastq
gzip sickle_cor_panda_B3.fastq
gzip sickle_cor_panda_C1.fastq
gzip sickle_cor_panda_C2.fastq
gzip sickle_cor_panda_C3.fastq
gzip sickle_cor_panda_D1.fastq
gzip sickle_cor_panda_D2.fastq
gzip sickle_cor_panda_D3.fastq
gzip sickle_cor_panda_E1.fastq
gzip sickle_cor_panda_E2.fastq
gzip sickle_cor_panda_E3.fastq
gzip sickle_cor_panda_F1.fastq
gzip sickle_cor_panda_F2.fastq
gzip sickle_cor_panda_F3.fastq

#############################################################################################

#### DAMe - using PCR replicates and read numbers to filter out bad reads

#1. place libraries in different folders
mkdir folder{A,B,C,D,E,F}{1,2,3}
mv sickle_cor_panda_A1.fastq.gz folderA1/
mv sickle_cor_panda_A2.fastq.gz folderA2/
mv sickle_cor_panda_A3.fastq.gz folderA3/
mv sickle_cor_panda_B1.fastq.gz folderB1/
mv sickle_cor_panda_B2.fastq.gz folderB2/
mv sickle_cor_panda_B3.fastq.gz folderB3/
mv sickle_cor_panda_C1.fastq.gz folderC1/
mv sickle_cor_panda_C2.fastq.gz folderC2/
mv sickle_cor_panda_C3.fastq.gz folderC3/
mv sickle_cor_panda_D1.fastq.gz folderD1/
mv sickle_cor_panda_D2.fastq.gz folderD2/
mv sickle_cor_panda_D3.fastq.gz folderD3/
mv sickle_cor_panda_E1.fastq.gz folderE1/
mv sickle_cor_panda_E2.fastq.gz folderE2/
mv sickle_cor_panda_E3.fastq.gz folderE3/
mv sickle_cor_panda_F1.fastq.gz folderF1/
mv sickle_cor_panda_F2.fastq.gz folderF2/
mv sickle_cor_panda_F3.fastq.gz folderF3/


#2. SORT

cd ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/folderA1
python /usr/local/bin/DAMe/bin/sort.py -fq sickle_cor_panda_A1.fastq.gz -p ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/Primers_COILeray.txt -t ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/Tags_300test_COIA.txt
cd ../folderA2
python /usr/local/bin/DAMe/bin/sort.py -fq sickle_cor_panda_A2.fastq.gz -p ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/Primers_COILeray.txt -t ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/Tags_300test_COIA.txt
cd ../folderA3
python /usr/local/bin/DAMe/bin/sort.py -fq sickle_cor_panda_A3.fastq.gz -p ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/Primers_COILeray.txt -t ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/Tags_300test_COIA.txt
cd ../folderB1
python /usr/local/bin/DAMe/bin/sort.py -fq sickle_cor_panda_B1.fastq.gz -p ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/Primers_COILeray.txt -t ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/Tags_300test_COIA.txt
cd ../folderB2
python /usr/local/bin/DAMe/bin/sort.py -fq sickle_cor_panda_B2.fastq.gz -p ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/Primers_COILeray.txt -t ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/Tags_300test_COIA.txt
cd ../folderB3
python /usr/local/bin/DAMe/bin/sort.py -fq sickle_cor_panda_B3.fastq.gz -p ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/Primers_COILeray.txt -t ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/Tags_300test_COIA.txt
cd ../folderC1
python /usr/local/bin/DAMe/bin/sort.py -fq sickle_cor_panda_C1.fastq.gz -p ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/Primers_COILeray.txt -t ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/Tags_300test_COIA.txt
cd ../folderC2
python /usr/local/bin/DAMe/bin/sort.py -fq sickle_cor_panda_C2.fastq.gz -p ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/Primers_COILeray.txt -t ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/Tags_300test_COIA.txt
cd ../folderC3
python /usr/local/bin/DAMe/bin/sort.py -fq sickle_cor_panda_C3.fastq.gz -p ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/Primers_COILeray.txt -t ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/Tags_300test_COIA.txt
cd ../folderD1
python /usr/local/bin/DAMe/bin/sort.py -fq sickle_cor_panda_D1.fastq.gz -p ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/Primers_COILeray.txt -t ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/Tags_300test_COIA.txt
cd ../folderD2
python /usr/local/bin/DAMe/bin/sort.py -fq sickle_cor_panda_D2.fastq.gz -p ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/Primers_COILeray.txt -t ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/Tags_300test_COIA.txt
cd ../folderD3
python /usr/local/bin/DAMe/bin/sort.py -fq sickle_cor_panda_D3.fastq.gz -p ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/Primers_COILeray.txt -t ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/Tags_300test_COIA.txt
cd ../folderE1
python /usr/local/bin/DAMe/bin/sort.py -fq sickle_cor_panda_E1.fastq.gz -p ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/Primers_COILeray.txt -t ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/Tags_300test_COIA.txt
cd ../folderE2
python /usr/local/bin/DAMe/bin/sort.py -fq sickle_cor_panda_E2.fastq.gz -p ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/Primers_COILeray.txt -t ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/Tags_300test_COIA.txt
cd ../folderE3
python /usr/local/bin/DAMe/bin/sort.py -fq sickle_cor_panda_E3.fastq.gz -p ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/Primers_COILeray.txt -t ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/Tags_300test_COIA.txt
cd ../folderF1
python /usr/local/bin/DAMe/bin/sort.py -fq sickle_cor_panda_F1.fastq.gz -p ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/Primers_COILeray.txt -t ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/Tags_300test_COIA.txt
cd ../folderF2
python /usr/local/bin/DAMe/bin/sort.py -fq sickle_cor_panda_F2.fastq.gz -p ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/Primers_COILeray.txt -t ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/Tags_300test_COIA.txt
cd ../folderF3
python /usr/local/bin/DAMe/bin/sort.py -fq sickle_cor_panda_F3.fastq.gz -p ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/Primers_COILeray.txt -t ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/Tags_300test_COIA.txt


# 3. Place PCR replicates in the same folder and rename to pool{1,2,3} for libraries A, B, C, D; pool{13,14,15} for library E;  pool{16,17,18} for library F
# The three PCR replicate folders (on which sort.py was run) have to be in the same folder (e.g. 'folderA') and named 'pool1', 'pool2', and 'pool3'. No other pool folders

cd ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/
mkdir folder{A,B,C,D,E,F}
mv folderA1 folderA/pool1
mv folderA2 folderA/pool2
mv folderA3 folderA/pool3
mv folderB1 folderB/pool1
mv folderB2 folderB/pool2
mv folderB3 folderB/pool3
mv folderC1 folderC/pool1
mv folderC2 folderC/pool2
mv folderC3 folderC/pool3
mv folderD1 folderD/pool1
mv folderD2 folderD/pool2
mv folderD3 folderD/pool3
mv folderE1 folderE/pool13
mv folderE2 folderE/pool14
mv folderE3 folderE/pool15
mv folderF1 folderF/pool16
mv folderF2 folderF/pool17
mv folderF3 folderF/pool18


# 4. Filter
# Filtering reads with low thresholds to get a feel for the data - 3 PCR replicates (pools), read present in min 2 pools, min 2 reads per pool, min 300 bp.  This step is slow (~ 45 mins per library).

cd ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/
python /usr/local/bin/DAMe/bin/filter.py -h

cd folderA
mkdir Filter_min2PCRs_min2copies_A
date; python /usr/local/bin/DAMe/bin/filter.py -psInfo ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/PSinfo_300test_COIA.txt -x 3 -y 2 -p 3 -t 2 -l 300 -o Filter_min2PCRs_min2copies_A; date

cd ../folderB
mkdir Filter_min2PCRs_min2copies_B
date; python /usr/local/bin/DAMe/bin/filter.py -psInfo ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/PSinfo_300test_COIA.txt -x 3 -y 2 -p 3 -t 2 -l 300 -o Filter_min2PCRs_min2copies_B; date

cd ../folderC
mkdir Filter_min2PCRs_min2copies_C
date; python /usr/local/bin/DAMe/bin/filter.py -psInfo ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/PSinfo_300test_COIA.txt -x 3 -y 2 -p 3 -t 2 -l 300 -o Filter_min2PCRs_min2copies_C; date

cd ../folderD
mkdir Filter_min2PCRs_min2copies_D
date; python /usr/local/bin/DAMe/bin/filter.py -psInfo ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/PSinfo_300test_COIA.txt -x 3 -y 2 -p 3 -t 2 -l 300 -o Filter_min2PCRs_min2copies_D; date

# Using different Psinfo files for libraries E and F.

cd ../folderE
mkdir Filter_min2PCRs_min2copies_E
date; python /usr/local/bin/DAMe/bin/filter.py -psInfo ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/PSinfo_300test_COIE.txt -x 3 -y 2 -p 3 -t 2 -l 300 -o Filter_min2PCRs_min2copies_E; date

cd ../folderF
mkdir Filter_min2PCRs_min2copies_F
date; python /usr/local/bin/DAMe/bin/filter.py -psInfo ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/PSinfo_300test_COIF.txt -x 3 -y 2 -p 3 -t 2 -l 300 -o Filter_min2PCRs_min2copies_F; date

## 5. Prepare the FilteredReads.fna file for Sumaclust clustering. changes the header lines on the fasta file

cd ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/
python /usr/local/bin/DAMe/bin/convertToUSearch.py -h

cd folderA/Filter_min2PCRs_min2copies_A
python /usr/local/bin/DAMe/bin/convertToUSearch.py -i FilteredReads.fna -lmin 300 -lmax 320
cd ../../folderB/Filter_min2PCRs_min2copies_B
python /usr/local/bin/DAMe/bin/convertToUSearch.py -i FilteredReads.fna -lmin 300 -lmax 320
cd ../../folderC/Filter_min2PCRs_min2copies_C
python /usr/local/bin/DAMe/bin/convertToUSearch.py -i FilteredReads.fna -lmin 300 -lmax 320
cd ../../folderD/Filter_min2PCRs_min2copies_D
python /usr/local/bin/DAMe/bin/convertToUSearch.py -i FilteredReads.fna -lmin 300 -lmax 320
cd ../../folderE/Filter_min2PCRs_min2copies_E
python /usr/local/bin/DAMe/bin/convertToUSearch.py -i FilteredReads.fna -lmin 300 -lmax 320
cd ../../folderF/Filter_min2PCRs_min2copies_F
python /usr/local/bin/DAMe/bin/convertToUSearch.py -i FilteredReads.fna -lmin 300 -lmax 320
# python /usr/local/bin/DAMe/bin/convertToUSearch.py -i FilteredReads.fna -lmin 300 -lmax 320 --usearch

## 6. Sumaclust clustering at 96% and 97% and convert Sumaclust output to table format

# sumaclust download from: https://git.metabarcoding.org/obitools/sumaclust/wikis/home
# wget https://git.metabarcoding.org/obitools/sumaclust/uploads/69f757c42f2cd45212c587e87c75a00f/sumaclust_v1.0.20.tar.gz
# tar -zxvf sumaclust_v1.0.20.tar.gz
# cd sumaclust_v1.0.20/
# make CC=clang # disables OpenMP, which isn't on macOS

cd ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/
~/src/sumaclust_v1.0.20/sumaclust -h
python /usr/local/bin/DAMe/bin/tabulateSumaclust.py -h

## 96% similarity
cd folderA/Filter_min2PCRs_min2copies_A/
~/src/sumaclust_v1.0.20/sumaclust -t 0.96 -e FilteredReads.forsumaclust.fna > OTUs_96_sumaclust.fna
python /usr/local/bin/DAMe/bin/tabulateSumaclust.py -i OTUs_96_sumaclust.fna -o table_300test_A_96.txt -blast

cd ../../folderB/Filter_min2PCRs_min2copies_B/
~/src/sumaclust_v1.0.20/sumaclust -t 0.96 -e FilteredReads.forsumaclust.fna > OTUs_96_sumaclust.fna
python /usr/local/bin/DAMe/bin/tabulateSumaclust.py -i OTUs_96_sumaclust.fna -o table_300test_B_96.txt -blast

cd ../../folderC/Filter_min2PCRs_min2copies_C/
~/src/sumaclust_v1.0.20/sumaclust -t 0.96 -e FilteredReads.forsumaclust.fna > OTUs_96_sumaclust.fna
python /usr/local/bin/DAMe/bin/tabulateSumaclust.py -i OTUs_96_sumaclust.fna -o table_300test_C_96.txt -blast

cd ../../folderD/Filter_min2PCRs_min2copies_D/
~/src/sumaclust_v1.0.20/sumaclust -t 0.96 -e FilteredReads.forsumaclust.fna > OTUs_96_sumaclust.fna
python /usr/local/bin/DAMe/bin/tabulateSumaclust.py -i OTUs_96_sumaclust.fna -o table_300test_D_96.txt -blast

cd ../../folderE/Filter_min2PCRs_min2copies_E/
~/src/sumaclust_v1.0.20/sumaclust -t 0.96 -e FilteredReads.forsumaclust.fna > OTUs_96_sumaclust.fna
python /usr/local/bin/DAMe/bin/tabulateSumaclust.py -i OTUs_96_sumaclust.fna -o table_300test_E_96.txt -blast

cd ../../folderF/Filter_min2PCRs_min2copies_F/
~/src/sumaclust_v1.0.20/sumaclust -t 0.96 -e FilteredReads.forsumaclust.fna > OTUs_96_sumaclust.fna
python /usr/local/bin/DAMe/bin/tabulateSumaclust.py -i OTUs_96_sumaclust.fna -o table_300test_F_96.txt -blast

# 97% similarity
cd ~/Xiaoyangmiseqdata/MiSeq_20170410/300Test/

cd folderA/Filter_min2PCRs_min2copies_A/
~/src/sumaclust_v1.0.20/sumaclust -t 0.97 -e FilteredReads.forsumaclust.fna > OTUs_97_sumaclust.fna
python /usr/local/bin/DAMe/bin/tabulateSumaclust.py -i OTUs_97_sumaclust.fna -o table_300test_A_97.txt -blast

cd ../../folderB/Filter_min2PCRs_min2copies_B/
~/src/sumaclust_v1.0.20/sumaclust -t 0.97 -e FilteredReads.forsumaclust.fna > OTUs_97_sumaclust.fna
python /usr/local/bin/DAMe/bin/tabulateSumaclust.py -i OTUs_97_sumaclust.fna -o table_300test_B_97.txt -blast

cd ../../folderC/Filter_min2PCRs_min2copies_C/
~/src/sumaclust_v1.0.20/sumaclust -t 0.97 -e FilteredReads.forsumaclust.fna > OTUs_97_sumaclust.fna
python /usr/local/bin/DAMe/bin/tabulateSumaclust.py -i OTUs_97_sumaclust.fna -o table_300test_C_97.txt -blast

cd ../../folderD/Filter_min2PCRs_min2copies_D/
~/src/sumaclust_v1.0.20/sumaclust -t 0.97 -e FilteredReads.forsumaclust.fna > OTUs_97_sumaclust.fna
python /usr/local/bin/DAMe/bin/tabulateSumaclust.py -i OTUs_97_sumaclust.fna -o table_300test_D_97.txt -blast

cd ../../folderE/Filter_min2PCRs_min2copies_E/
~/src/sumaclust_v1.0.20/sumaclust -t 0.97 -e FilteredReads.forsumaclust.fna > OTUs_97_sumaclust.fna
python /usr/local/bin/DAMe/bin/tabulateSumaclust.py -i OTUs_97_sumaclust.fna -o table_300test_E_97.txt -blast

cd ../../folderF/Filter_min2PCRs_min2copies_F/
~/src/sumaclust_v1.0.20/sumaclust -t 0.97 -e FilteredReads.forsumaclust.fna > OTUs_97_sumaclust.fna
python /usr/local/bin/DAMe/bin/tabulateSumaclust.py -i OTUs_97_sumaclust.fna -o table_300test_F_97.txt -blast

#### ***** Using CROP makes 346 97% OTUs, versus 391 sumaclust 97% OTUs
cd folderA/Filter_min2PCRs_min2copies_A/
~/src/CROP/CROP -i FilteredReads.forsumaclust.fna -o OTUs_CROP_97.out -s -e 2000 -b 185 -z 300 -r 0
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

# START HERE

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
