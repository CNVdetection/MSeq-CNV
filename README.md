# MSeq-CNV
MSeq-CNV: accurate detection of Copy Number Variation from Multiple Sequencing samples

By: 
Seyed Amir Malekpour, Hamid Pezeshk and Mehdi Sadeghi

Contact: a.malekpour@ut.ac.ir, pezeshk@khayam.ut.ac.ir, pezeshk@ut.ac.ir, sadeghi@nigeb.ac.ir

MSeq-CNV can be applied for detecting the recurrent genome-wide CNVs from NGS data in the diploid genome of human and other organisms, as well. The input NGS data for the MSeq-CNV are possibly the mate pair reads which are collected from sequencing with multiple platforms, multiple individuals and experimental conditions.
Run MAIN.m MATLAB script, to call genome-wide deletions and duplications using MSeq-CNV. However, prior to running this script, the input data and a few variables have to be specified in the MATLAB workspace:

-	Input data: the input data should be saved in a MATLAB matrix, named “Info”. Each row of this matrix corresponds to a mate pair, with the four entries:
1-	Leftmost position of where the first read maps to the reference, 
2-	Leftmost position of where the next read maps to the reference, 
3-	Mate pair insertion size, after mapping the mate pair to the reference genome
4-	Sample number (mate pairs which are generated from the same sample genome, receive an identical number). For example, if six sample genomes are analyzed using MSeq-CNV, samples are given numbers 1, 2, 3, 4, 5 and 6, in order. 
Items 1, 2 and 3 are extracted from SAM or BAM files and are then used in MSeq-CNV.


There are other variables which are given values in the first line of the MAIN script: 
-	Chromosome size (in bp)
-	Sample Size (with a default value of 6)
-	Segment Size (with a default value of 110 bp)
-	Initial estimate of Lambda (The average number of reads which are expected to map to a genomic segment with diploid state, the default value of this variable is 42)
-	T (this variable identifies the number of genomic segments which are analyzed using MSeq-CNV, with a default value of 250,000)
-	Insertion Size (This is the average mate pair insertion size in the clone library, with a default value of 150 bp)
-	Read Length (The length of each read in a mate pair, with a default value of 35 bp)
-	A sample of the observed mate pair insertion sizes in the clone library (required for estimating the empirical distribution of the mate pair insertion sizes in the clone library) 
-	The expected minimum of the mate pair insertion sizes in the clone library
-	The expected maximum of the mate pair insertion sizes in the clone library 


After running MAIN.m, genome-wide copy gain and copy loss regions are saved in a MATLAB cell i.e. “CNVs”. For example, genomic regions with copy number variation in the first sample genome are saved in MATLAB variable CNVs{1}. Each row of this matrix represents a different CNV call with its start position, end position and copy number estimation, in order.
