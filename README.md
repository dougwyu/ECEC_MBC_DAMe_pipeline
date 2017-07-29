# ECEC_MBC_DAMe_pipeline

AdapterRemoval_sickle_BayesHammer_PandaSeq_DAMe

See Wiki page for assumed folder structure

This pipeline is tested for macOS and processes Illumina HiSeq/MiSeq files for metabarcoding. Software installation information is in the shell script as comments.

Create merged reads (following Schirmer et al. 2015, NAR)
1)  AdapterRemoval to remove adapters
2)  sickle for trimming
3)  BayesHammer in spades.py (--error-correction-only) for denoising
4)  PandaSeq for read merging

Filter and Diagnosis the merged reads using the DAMe pipeline (Zepeda-Mendoza et al. 2016, BMC Research Notes)
1)  sort.py # remove tags from 

Note that the DAMe pipeline assumes that each set of samples (a 96-well plate) has been separately PCRd three times.  DAME keeps reads that appear in ≥M PCRs (e.g. ≥2 of 3 total PCRs) and in each PCR appears in ≥N copies (e.g. ≥4 copies). DAMe also works best if both the forward and reverse reads have been tagged, using the same tag (what we call 'twin tags').
