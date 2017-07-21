#!/bin/bash
set -e
set -u
set -o pipefail
#######################################################################################
#######################################################################################
# a shell script to loop through a set of fastq files and assemble all the mitogenomes
#######################################################################################
#######################################################################################

# Usage: bash _loop_mitogenome-assembly.sh MINCONTIG
# e.g. bash _loop_mitogenome-assembly.sh 6000
# try this:  bash _loop_mitogenome-assembly.sh 6000 | tee _loop_mitogenome-assembly$(date).log

# Interactive procedure:  run these commands before running _loop_mitogenome-assembly.sh interactively
# interactive -x -R "rusage[mem=20000]" -M 22000  # alternative is:  interactive -q interactive-lm
# module load gcc/4.9.3; module load java/jre/1.6.0-24; module load blast/ncbi-2.2.31+; module load samtools/1.2; module load bwa/0.7.9a
# upload _loop_mitogenome-assembly.sh in Project_TOG (e.g. ~/greenland_2016/Copenhagen_mtgenome_refs/20160610_hiseq6b/Project_TOG/)
# cd ~/greenland_2016/Copenhagen_mtgenome_refs/20160610_hiseq6b/Project_TOG

# Batch procedure:  bsub script and submit to long-lm queue
	#######################################################################################
	#######################################################################################
	#!/bin/sh
	#BSUB -q long-lm     #debug, short (24 hours), medium (24 hours), long (168 hours = 7 days)
	#BSUB -J mitoloop
	#BSUB -oo mitoloop.out
	#BSUB -eo mitoloop.err
	#BSUB -R "rusage[mem=20000]"
	#BSUB -M 22000
	#BSUB -B        #sends email to me when job starts
	#BSUB -N        # sends email to me when job finishes
	#
	# . /etc/profile
	# module purge
	# module load gcc/4.9.3
	# module load java/jre/1.6.0-24
	# module load blast/ncbi-2.2.31+
	# module load samtools/1.2
	# module load bwa/0.7.9a
	#
	# bash _loop_mitogenome_assembly_20160708.sh 6000
	#######################################################################################
	#######################################################################################


PIPESTART=$(date)

# upload _loop_mitogenome_assembly.sh to folder that contains all the genome skims (Project_TOG)
HOMEFOLDER=$(pwd)

echo "Home folder is "${HOMEFOLDER}""

# set variables
MINCONTIG=$1 # min contig length after assembly
EVALUE="1e-7"
SMALLCONTIGS="3000"
INDEX=1
if [ ! -e species_without_mitogenome_contigs.txt ] # if species_without_mitogenome_contigs.txt does not exist
then
	touch species_without_mitogenome_contigs.txt
fi

if [ ! -d pilon_outputs ] # if directory pilon_outputs does not exist
then
	mkdir pilon_outputs
fi

# # read in bash array of folders from a list that i make by hand
	# sample_info=104spp_fastq_files_array.txt # folder names in column 2, 5 unique samples
	# sample_names=($(cut -f 2 "$sample_info" | uniq))
	# # echo ${sample_names[@]} # echo all array elements
	# echo "There are" ${#sample_names[@]} "folders that will be processed." # echo number of elements in the array

# read in folder list and make a bash array
find * -maxdepth 0 -type d -name "Sample*" > folderlist.txt  # find all folders starting with "Sample*"
sample_info=folderlist.txt # put folderlist.txt into variable
sample_names=($(cut -f 1 "$sample_info" | uniq)) # convert variable to array this way
# echo ${sample_names[@]} # echo all array elements
echo "There are" ${#sample_names[@]} "folders that will be processed." # echo number of elements in the array

for sample in ${sample_names[@]}  # ${sample_names[@]} is the full bash array
do
	cd "${HOMEFOLDER}"
	echo "Now on species" ${INDEX} of ${#sample_names[@]}". Moved back to starting directory:"
	INDEX=$((INDEX+1))
	pwd
	FOLDER="${sample}" # this sets the folder name to the value in the bash array (which is a list of the folders)
	# SPECIES="$(find "${sample}" -name "*R1_001.fastq.gz" | xargs basename -s "_R1_001.fastq.gz")"  # -s used in macOS version
		# the output of the find command is piped to basename and basename strips the path and the suffix, and the stripped name is saved to $SPECIES
		# the output of the $(find | xargs basename) is wrapped with $() for command substitution, and wrapped in "" for robustness
		# the name is saved into $SPECIES
	SPECIES0=$(find "${sample}" -name "*R1_001.fastq.gz") # the piped xargs version returns the suffix on GRACE, not the filename.  weird
	SPECIES=$(basename ${SPECIES0} _R1_001.fastq.gz) # so i have split into two commands, which works on the command line
	SPECIES0= # clear the variable SPECIES0
		# this is complicated because the original version using xargs doesn't work on GRACE
		# the find command finds the *R1_001.fastq.gz filename inside folder $sample and saves it to $SPECIES0
		# basename strips the path and the suffix _R1_001.fastq.gz and saves to $SPECIES (i used to do this by hand)

	echo "**** Minimum acceptable contig length is" $MINCONTIG "bp"
	echo "**** blast evalue to filter-in putative mitogenomic contigs" is $EVALUE
	echo "**** Working folder is" $FOLDER
	echo "**** Species name is" $SPECIES


	mkdir ${SPECIES}_working
	# echo "**** copying fastq.gz files to working folder"
	# cp "${FOLDER}/${SPECIES}_R1_001.fastq.gz" "${SPECIES}_working/"
	# cp "${FOLDER}/${SPECIES}_R2_001.fastq.gz" "${SPECIES}_working/"
	# echo "**** fastq.gz files moved to working folder"
	echo "moving to ${SPECIES}_working"
	cd ${SPECIES}_working; pwd

	#### denoise using bfc
	echo "**** start of bfc denoising"
	mkdir bfc
	~/scripts/bfc/bfc -s 3g -t16 "${FOLDER}/${SPECIES}_R1_001.fastq.gz" | gzip -1 > bfc/${SPECIES}_R1_bfc.fastq.gz
	~/scripts/bfc/bfc -s 3g -t16 "${FOLDER}/${SPECIES}_R2_001.fastq.gz" | gzip -1 > bfc/${SPECIES}_R2_bfc.fastq.gz
	echo "**** end of bfc denoising"
	#### use Trimmomatic to remove illumina adapters
	# use TruSeq3-PE-2.fa file, which is the latest version (http://seqanswers.com/forums/archive/index.php/t-50419.html)
	echo "**** start of Trimmomatic trimming"
	cd bfc/ ; pwd
	java -jar ~/scripts/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 16 -phred33 -basein ${SPECIES}_R1_bfc.fastq.gz -baseout ${SPECIES}_bfc_trim.fastq.gz ILLUMINACLIP:~/scripts/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10:1:true LEADING:3 TRAILING:3 MINLEN:110 &> myjob.log
	echo "**** end of Trimmomatic trimming"

	#### move trimmed fastq files to Trimmed/ folder
	cd ../ ; pwd
	mkdir Trimmed
	mv bfc/*P.fastq.gz Trimmed/
	cd Trimmed/ ; pwd

	# echo "**** start of seqtk merge interleaving"
	~/scripts/seqtk/seqtk mergepe ${SPECIES}_bfc_trim_1P.fastq.gz ${SPECIES}_bfc_trim_2P.fastq.gz > ${SPECIES}_bfc_trim_12P_merged.fastq.gz
	# echo "**** end of seqtk mergepe interleaving"

	#### move the new, merged fasta files into Assembly/ and run
	################## SPAdes assembly ##################################
	# Nurk, S., Bankevich, A., Antipov, D., Gurevich, A., Korobeynikov, A., Lapidus, A., Prjibelsky, A., Pyshkin, A., Sirotkin, A., Sirotkin, Y., Stepanauskas, R., McLean, J., Lasken, R., Clingenpeel, S.R., Woyke, T., Tesler, G., Alekseyev, M.A. & Pevzner, P.A. (2013) Assembling Genomes and Mini-metagenomes from Highly Chimeric Reads. Research in Computational Molecular Biology Lecture Notes in Computer Science. pp. 158–170. Springer Berlin Heidelberg, Berlin, Heidelberg.
	# http://lab.loman.net/page7/
	echo "**** start of assembly"
	cd ../ # move out of Trimmed

	fn1="../Trimmed/${SPECIES}_bfc_trim_1P.fastq.gz"
	fn2="../Trimmed/${SPECIES}_bfc_trim_2P.fastq.gz"
	fn12="../Trimmed/${SPECIES}_bfc_trim_12P_merged.fastq.gz"
	out="Assembly"

	~/scripts/SPAdes-3.8.2-Linux/bin/spades.py --only-assembler --memory 46 -k 21,33,55 --pe1-1 ${fn1} --pe1-2 ${fn2} -o ${out}
	# ~/scripts/SPAdes-3.8.2-Linux/bin/spades.py --careful --memory 46 -k 21,33,55,77 --pe1-1 ${fn1} --pe1-2 ${fn2} -o ${out}  #
	# ~/scripts/SPAdes-3.8.2-Linux/bin/spades.py --only-assembler --memory 46 -k 21,33,55,77 --pe1-12 ${fn12} -o ${out}
	# --careful (= MismatchCorrector) adds 12 hours to the pipeline, and pilon does the same thing almost immediately

	# --meta # for metagenomic data
	# --plasmid # for plasmid assembly
	# --careful # runs MismatchCorrector, can't be used with metagenomic data
	# --continue # continues a crashed assembly from last kmer (k) completed
	# --only-assembler # runs assembler only, no error-correction
	# --only-error-correction # runs error-correction module only
	# --pe1-1 ${fn1} # forward reads, first library
	# --pe1-2 ${fn2} # reverse reads, first library
	# -m 22 # 22 GB max memory allowed to be used
	# -k 21,33,55,77 # kmer steps
	# -k 21,33,55,77,99,127 # kmer options for 250 PE libraries

	NOW=$(date)
	echo "**** ended SPAdes assembly at $NOW"
	#### filter through the contigs to find those that are mitochondrial
	echo "**** removing contigs < ${SMALLCONTIGS} bp and blastn filtering"
	# remove small contigs (< $SMALLCONTIGS bp)
	~/scripts/seqtk/seqtk seq -L ${SMALLCONTIGS} "${out}"/contigs.fasta > contig_min"${SMALLCONTIGS}".fasta

	blastn -db ~/scripts/mitoarthropoda/mitoarthropoda_20160625.fasta -query "${out}"/contig_min"${SMALLCONTIGS}".fasta -out ${SPECIES}_SPAdescontigs.blastn -num_threads 16 -task blastn -evalue $EVALUE -max_target_seqs 1 -outfmt 6 -dust yes
	# default value 1e-7


	########################################################################################################################################################################
	# by hand
	# interactive
	# cd greenland_2016/Copenhagen_mtgenome_refs/20160610_hiseq6b/Project_TOG; module load blast/ncbi-2.2.31+
	# cd to correct Working folder
	# ~/scripts/seqtk/seqtk seq -L 7000 contigs.fasta > contig_min7000.fasta; blastn -db ~/scripts/mitoarthropoda/mitoarthropoda_20160625.fasta -query contig_min7000.fasta -out SPAdescontigs.blastn -num_threads 16 -task blastn -evalue 1e-7 -max_target_seqs 1 -outfmt 6 -dust yes; cut -f 1 SPAdescontigs.blastn | sort -u > SPAdescontigs.blastn.list; ~/scripts/seqtk/seqtk subseq contigs.fasta SPAdescontigs.blastn.list > SPAdescontigs_min7000.fasta; ~/scripts/seqtk/seqtk seq -L 7000 scaffolds.fasta > scaffolds_min7000.fasta; blastn -db ~/scripts/mitoarthropoda/mitoarthropoda_20160625.fasta -query scaffolds_min7000.fasta -out SPAdesscaffolds.blastn -num_threads 16 -task blastn -evalue 1e-7 -max_target_seqs 1 -outfmt 6 -dust yes; cut -f 1 SPAdesscaffolds.blastn | sort -u > SPAdesscaffolds.blastn.list; ~/scripts/seqtk/seqtk subseq scaffolds.fasta SPAdesscaffolds.blastn.list > SPAdesscaffolds_min7000.fasta
	########################################################################################################################################################################

	########################################################################################################################################################################
	#  check if ${SPECIES}_IDBAcontigs.blastn is empty. if empty, use continue to jump out of this iteraction of the loop
	BLASTOUTPUT="${SPECIES}_SPAdescontigs.blastn"
	if [ ! -s $BLASTOUTPUT ]  # if $BLASTOUTPUT does not have any content (file size is 0)
	then
	  cd "${HOMEFOLDER}"
		echo "${FOLDER} produced no mitogenome contigs" >> species_without_mitogenome_contigs.txt  # append notice to file
		continue # jump out of this interation of the loop
	fi
	########################################################################################################################################################################

# keep only first column, sort by keep only uniques
	cut -f 1 ${SPECIES}_SPAdescontigs.blastn | sort -u > ${SPECIES}_SPAdescontigs.blastn.list

	# filter in mitochondrial reads to a new fasta file (should use seqtk subseq)
	# look through these contigs for mitogenomes (MITOS webserver)
	~/scripts/seqtk/seqtk subseq "${out}"/contigs.fasta ${SPECIES}_SPAdescontigs.blastn.list > ${SPECIES}_SPAdescontigs_min${MINCONTIG}.fasta
	# output file is stored in Trimmed, not spadesfiles!!
	echo "**** end of blastn filtering"

	########################################################################################################################################################################
	#  check if ${SPECIES}_SPAdescontigs_min${MINCONTIG}.fasta is empty. If empty, use continue to jump to next species
	SEQTKOUTPUT="${SPECIES}_SPAdescontigs_min${MINCONTIG}.fasta"
	if [ ! -s $SEQTKOUTPUT ]  # if $BLASTOUTPUT does not have any content (file size is 0)
	then
		cd "${HOMEFOLDER}"
		echo "${FOLDER} produced no mitogenome contigs" >> species_without_mitogenome_contigs.txt  # append notice to file
		continue # jump out of this interation of the loop
	fi
	########################################################################################################################################################################

	#### bwa to samtools to Pilon for consensus contig.fa ####
	# index the reference genome
	echo "**** start of bwa indexing of contig.fa"
	echo "**** If pipeline crashes here, could be due to no sequences in: ${SPECIES}_SPAdescontigs_min${MINCONTIG}.fasta"
	# uses bwa 0.7.15-r1142-dirty, installed locally on 1 Jul 2016 from github
	~/scripts/bwa-gh/bwa index -a is ${SPECIES}_SPAdescontigs_min${MINCONTIG}.fasta  # <-a is>   uses the is format, which is appropriate for short genomes, like mitochondria

	# uses bwa 0.7.12-r1039, installed locally on 1 Jul 2016 from sourceforge
	# ~/scripts/bwa-sf/bwa index -a is ${SPECIES}_IDBAcontigs_min${MINCONTIG}.fasta  # <-a is>   uses the is format, which is appropriate for short genomes, like mitochondria
	# ~/scripts/bwa-sf/bwa index -a is scaffold.fa # idba produces scaffolds, and i want to save these for possible minHash approach
	echo "**** end of bwa indexing of contig.fa"
	# align the reads and save input to a sam file
	# I am using the post-bfc-denoised, post-trimmomatic-trimmed sequences, but before the usearch -maxeee 1.0 filtering
	# Of course, i could use the post usearch -maxee filtered read files, but i thought i would use the bigger dataset

	echo "**** start of bwa mem mapping, and sam to bam conversion of contig.fa"
	# uses bwa 0.7.15-r1142-dirty, installed locally on 1 Jul 2016 from github, piped to samtools without making samfile
	~/scripts/bwa-gh/bwa mem -t 16 -v 3 ${SPECIES}_SPAdescontigs_min${MINCONTIG}.fasta ../${SPECIES}_bfc_trim_1P.fastq.gz ../${SPECIES}_bfc_trim_2P.fastq.gz | samtools view -b -S -o ${SPECIES}_bfc_trim_P.bam -
		# -b: indicates that the output is BAM
		# -S: indicates that the input is SAM
		# -o: specifies the name of the output file
	echo "**** end of bwa mem mapping, and sam to bam conversion of contig.fa"

	# uses bwa 0.7.12-r1039, installed locally on 1 Jul 2016 from sourceforge
	# ~/scripts/bwa-sf/bwa mem -t 16 -v 3 ${SPECIES}_SPAdescontigs_min${MINCONTIG}.fasta ../Trimmed/${SPECIES}_bfc_trim_1P.fastq.gz ../Trimmed/${SPECIES}_bfc_trim_2P.fastq.gz | samtools view -b -S -o ${SPECIES}_bfc_trim_P.bam

	# samtools view -b -S -o ${SPECIES}_bfc_trim_P.bam ${SPECIES}_bfc_trim_P.sam
	# echo "**** end of samtools sam to bam conversion"

	echo "**** start of samtools bam sorting of contig.fa"
	samtools sort ${SPECIES}_bfc_trim_P.bam ${SPECIES}_bfc_trim_P.sorted
	echo "**** end of samtools bam sorting of contig.fa"

	echo "**** start of samtools bam indexing of contig.fa"
	samtools index ${SPECIES}_bfc_trim_P.sorted.bam  # makes a *.bai file
	echo "**** end of samtools bam indexing of contig.fa"

	# Pilon correction
	echo "**** start of pilon correction of contig.fa"
	# -Xmx24G requests 24GB memory, can add <--vcf> if
	java -Xmx24G -jar ~/scripts/pilon/pilon-1.18.jar --genome ${SPECIES}_SPAdescontigs_min${MINCONTIG}.fasta --frags ${SPECIES}_bfc_trim_P.sorted.bam --outdir pilon --changes &> pilon_contig.log
	mv pilon_contig.log pilon/pilon_contig.log

	mv pilon/pilon.fasta pilon/${SPECIES}_contig_pilon.fasta
	mv pilon/pilon.changes pilon/${SPECIES}_contig_changes_pilon.txt
	mv pilon/pilon_contig.log pilon/${SPECIES}_contig_pilon.log

	cp pilon/${SPECIES}_contig_* "${HOMEFOLDER}"/pilon_outputs

	echo "**** ${SPECIES}_contig_pilon.fasta"
	echo "**** ${SPECIES}_contig_changes_pilon.txt"
	echo "**** ${SPECIES}_contig_pilon.log"

	echo "**** end of pilon assembly correction of contig.fa"

	# format of ${SPECIES}_pilon_changes.txt:  <Original Scaffold Coordinate> <New Scaffold Coordinate> <Original Sequence> <New Sequence>
	# contig-100_0:16312 contig-100_0_pilon:16312 A G
	# contig-100_1:180-183 contig-100_1_pilon:180 AAAG .
	# contig-100_1:1927 contig-100_1_pilon:1923-1924 . GA
	# contig-100_1:9624 contig-100_1_pilon:9622 T C


########################################################################################################################################################################
#  if spades produced a scaffold.fa, run bwa and pilon on it for a possible minHash approach
	if [ -e spadesfiles/scaffolds.fasta ]
	then
		gzip spadesfiles/scaffolds.fasta
		cp spadesfiles/scaffolds.fasta.gz "${HOMEFOLDER}"/pilon_outputs/"${SPECIES}"_scaffold.fa.gz
	fi
	#################################################################################################
	# send pilon-corrected fasta(s) to the MITOS webserver and download the GFF file to add the annotation
	# save the annotated fasta file in Geneious

### activate when i am happy with the loop code
	cd "${HOMEFOLDER}"
	rm -rf "${SPECIES}"_working # remove the working directory to make space

done

echo "Pipeline started at $PIPESTART"
NOW=$(date)
echo "Pipeline ended at   $NOW"
echo "**** end of assembly pipeline. finished assemblies are in 20160610_hiseq6b/Project_TOG/pilon_outputs folder"
