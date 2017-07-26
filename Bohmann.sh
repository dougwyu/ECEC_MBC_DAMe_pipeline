P=9
for i in `seq 1 $P`
do
cd /groups/hologenomics/kbohmann/data/Vampire/16s/pool${i}
python /groups/hologenomics/software/DAMe/v0.9/bin/sort.py -fq Pool${i}_merged.fastq -p /groups/hologenomics/kbohmann/data/Vampire/16s/Primers_16sMam.txt -t /groups/hologenomics/kbohmann/data/Vampire/16s/Tags_16sMam.txt
done

#sorting the summary counts file (it’s on my list to get it incorporated in the sort.py - it’s nicer when it’s sorted :) )

P=9
for i in `seq 1 $P`
do
cd /groups/hologenomics/kbohmann/data/Vampire/16s/pool${i}
head -1 SummaryCounts.txt > SummaryCounts_sorted.txt
tail -n +2 SummaryCounts.txt | sed "s/Tag//g" | sort -k1,1n -k2,2n | awk 'BEGIN{OFS="\t";}{$1="Tag"$1;$2="Tag"$2; print $0;}' >> SummaryCounts_sorted.txt
cd ..
done

#making the tag combinations overview
P=9
for pool in `seq 1 $P`; do
python /groups/hologenomics/software/DAMe/v0.9/bin/splitSummaryByPSInfo.py -p PSinfo_vampire_16s.txt -l $pool -s pool$pool/SummaryCounts_sorted.txt -o pool$pool/SummaryCounts_split.txt
done

#filter
python /groups/hologenomics/software/DAMe/v0.9/bin/filter.py -psInfo PSinfo_vampire_16s.txt -x 3 -y 2 -p 9 -t 20 -l 80

python /groups/hologenomics/software/DAMe/v0.9/bin/plotLengthFreqMetrics_perSample.py -f FilteredReads.fna -p /groups/hologenomics/kbohmann/data/Vampire/16s/PSinfo_vampire_16s.txt -n 9

#assessing clustering parameters
xsbatch -c 4 --mem-per-cpu=5000 --time=6:00:00 -- python /groups/hologenomics/software/DAMe/v0.9/bin/assessClusteringParameters.py -i FilteredReads.fna -mint 0.8 -minR 0.6 -step 0.05 -t 4 -o 16sclusterassess_mint08_minR06_step005.pdf

python /groups/hologenomics/software/DAMe/v0.9/bin/convertToUSearch.py -i FilteredReads.fna -lmin 80 -lmax 100

sumaclust -t 0.98 -e FilteredReads.forsumaclust.fna > OTU98/OTU98_vamps_16s.fna

#make OTUxSample table
python /groups/hologenomics/software/DAMe/v0.9/bin/tabulateSumaclust.py -s 50000 -i OTU98_vamps_16s.fna -o supertable_Vamps_16s_98.txt -blast
