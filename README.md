# ECEC_MBC_DAMe_pipeline

See Wiki page for assumed folder structure.  

**This pipeline is designed for my laboratory to use, so clone and employ in your environment with care!**

This pipeline is tested for macOS and processes Illumina HiSeq/MiSeq files for metabarcoding. Software installation information is in the shell script as comments.

Create merged reads (following Schirmer et al. 2015, NAR)
1)  AdapterRemoval to remove adapters
2)  sickle for trimming
3)  BayesHammer in spades.py (--error-correction-only) for denoising
4)  PandaSeq for read merging

Filter, Diagnose, and Cluster the merged reads using the DAMe pipeline (Zepeda-Mendoza et al. 2016, BMC Research Notes)
1)  sort.py # remove tags from merged reads and count them up
2)  some diagnostics to detect incorrect sample tagging and PCR problems
3)  filter.py # keep reads that exceed Minimum_PCR and Minimum_copy_number_per_PCR thresholds, defined by the user
4)  Renkonen Similarity Index analysis to check for PCR errors
5)  sumatra analysis to inform OTU similarity-threshold choice
6)  vsearch -uchime_denovo chimera removal
7)  sumaclust OTU clustering and creation of the OTU table
8)  RDP Classifier taxonomic assignment
9)  Filter OTU tables and OTU representative sequences fasta by taxonomy (only Arthropoda ≥ 0.80 probability)

Prepare OTU tables for  
1) Merge OTU table and taxonomies
2) Filter likely contaminant OTUs as indicated by negative and positive controls
3) Filter out artefactually small OTUs via phyloseq method (https://github.com/joey711/phyloseq)
4) Create and save full OTU tables (OTU x samples) and community-analysis tables (sample X OTUs)
5) Simple NMDS ordination

Note that the DAMe pipeline assumes that each set of samples (a 96-well plate) has been separately PCRd three times.  DAME keeps reads that appear in ≥M PCRs (e.g. ≥2 of 3 total PCRs) and in each PCR appears in ≥N copies (e.g. ≥4 copies). DAMe also works best if both the forward and reverse reads have been tagged, using the same tag (what we call 'twin tags').
