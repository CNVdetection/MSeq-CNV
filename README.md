# MSeq-CNV
MSeq-CNV: accurate detection of Copy Number Variation from Sequencing of Multiple samples

By: 
Seyed Amir Malekpour, Hamid Pezeshk and Mehdi Sadeghi

Contact: a.malekpour@ut.ac.ir, pezeshk@khayam.ut.ac.ir, pezeshk@ut.ac.ir, sadeghi@nigeb.ac.ir

MSeq-CNV can be applied for detecting recurrent genome-wide CNVs from NGS data in the diploid genome of human and other organisms, as well. 
The input NGS data for the MSeq-CNV are possibly the mate pair reads which are collected from sequencing with multiple platforms, multiple individuals and experimental conditions.
MSeq-CNV programs are written in R language. Run MseqCNV.R script, to call genome-wide deletions and duplications using MSeq-CNV. 
However, prior to running this script, the input data and a few variables have to be specified in the beginning of this script:


Input data: the input data should be saved in a .csv file. Each row of this matrix corresponds to a mate pair, with the four entries:

1- Leftmost position of where the first read maps to the reference,
2- Mate pair insertion size, after mapping the mate pair to the reference genome
3- Mate pair flag which indicates its mapping properties
4- Mate pair mapping quality
5- Read number i.e. 1st or 2nd in a mate pair 
6- Sample number (mate pairs which are generated from the same sample genome, receive an identical number). 
For example, if six sample genomes are analyzed using MSeq-CNV, samples are given numbers 1, 2, 3, 4, 5 and 6, in order.

Items 1, 2, 3, 4 and 5 are extracted from SAM or BAM files and are then used in MSeq-CNV.


To detect CNVs using MseqCNV, first download the constructed sample genome data from our page i.e. SampleData.csv. 

There are variables which are given values in the first line of the MseqCNV.R script: 
- Path to the .csv data file
- Chromosome size (in bp)
- Sample Size (with a default value of 40)
- Segment Size (with a default value of 150 bp)
- Initial estimation of Lambda (The average number of reads which are expected to map to a genomic segment with diploid state, the default value of this variable is 13)
- T (this variable identifies the number of genomic segments which are analyzed using MSeq-CNV, with a default value of 1,000)
- Insertion Size (This is the average mate pair insertion size in the clone library, with a default value of 200 bp)
- Read Length (The length of each read in a mate pair, with a default value of 50 bp)
- Number of iterations in EM algorithm

e.g. 
Pos_Size <- read.csv("D:\\SampleData.csv")

chromosome_size =5000000

sample_size=40

segment_size=150

lambda_initial_estimate <-  13

T=1000

Insertion_Size=200

Read_Length=50

iteration=3



After running MseqCNV.R, genome-wide copy gain and copy loss regions are saved in a matrix named CNVs.
Each row of this matrix represents a different CNV call with its sample number, start position, end position and copy number estimation.
