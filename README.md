transcriptome
=============

Pipeline for Illumina transcriptome data.
Sonal Sinal & Ke Bi
July 31 2013


**I am still working on writing it. If you have any questions please contact me (Ke) directly**

If you use this pipeline please cite: Singhal 2013 De novo transcriptomic analyses for non-model organisms: an evaluation of methods across a multi-species data set. Molecular Ecology Resources DOI: 10.1111/1755-0998.12077

****************************************************************
1-pre-cleanup.pl: format raw sequence reads for 2-scrubReads.pl

The input is the raw sequence libraries directly from Illumina which should be something like:
KB1_ATCACG_L002_R1_001.fasta.gz

If the input is something like above, then the output (in a folder called "pre-clean") are:
KB1_R1.fq


*************************************************************** 
2-scrubReads.pl: remove adapters/contaminations/low complexities and merge overlapping reads.

The input is the output from 1-pre-cleanup.pl 
for example: KB1_R1.fq
The output contains six txt files:
KB1_1_final.txt (cleaned left reads)
KB1_2_final.txt (cleaned right reads)
KB1_u_final.txt (cleaned, merged and unpaired reads)
KB1.contam.out (reads derived from contamination)
KB1.duplicates.out (duplicated reads)
KB1.lowComplexity.out (low complex reads)

Examples for library info and adapter files can be found in the folder named "Sample" in the same repository

Note for adapter file: if you have only P7 index adapters just leave the P5 column "NA" or blank for every library, or just exclude the P5 column.  
******************************************************************

3-generateAssemblies_LOCAL.pl: running ABySS (exon capture data only) and Trinity (cDNA data only) using local computers/clusters

ABySS uses OpenMP for parallelization and requires the Boost C++ libraries. It should be also built using the sparsehash library to reduce memory usage.
For details please visit http://seqanswers.com/wiki/ABySS

******************************************************************

3-generateAssemblies_XSEDE.pl: a shell script for Trinity job submission on the Blacklight.

******************************************************************




