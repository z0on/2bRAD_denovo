Whole genome de novo genotyping with 2bRAD
------------------------------------------

Mikhail Matz, matz@utexas.edu

2bRAD has been described in Wang et al 2012 
http://www.nature.com/nmeth/journal/v9/n8/abs/nmeth.2023.html 

2bRAD features a very simple library prep protocol with no intermediate purification stages. It involves triple barcode scheme, two being standard Illumina barcodes and one in-read ligated barcode for pooling libraries in 12-plexes midway through library prep to further minimize prep costs. 2bRAD generates 36-base tags that can be sequenced on either strand. In a full-representation version it results in one high-quality genotyped tag every 2.5-3.5 kb on average. With reduced-representation with modified ligated adapters, the number of genotyped tags can be reduced 4- or 16-fold, to save sequencing costs in applications not requiring dense genome sampling (such as population structure analysis or linkage mapping). 

Unique features of 2bRAD bioinformatics are removal of PCR duplicates based on degenerate tag (allowing for 128-fold dynamic range of sequencing depth per allele), use of direct-reverse strand ratio as a quality fltering parameter, as well filtering and quality assessment based on genotyping replicates.

de novo 2bRAD pipeline is similar in ideology to STACKS, but takes advantage of the fact that 2bRAD sequences both strands. The genotype caller (haplocall_denovo.pl) emits VCF 4.1 formatted files (both for whole tags and individual SNPs). This allows for easy interface with downstream analysis software, most importantly vcftools that can be used for calculating Fst, nucleotide diversity, Tajima's D, heterozygosity, and linkage disequilibrium. 

This repository contains the lab protocol for sample preparation, as well as scripts and walkthrough (2bRAD_denovo_README.txt) for:
- trimming and quality filtering;
- removing PCR duplicates;
- assembling loci;
- calling variants (SNP-wise and haplotype-wise);
- recalibrating quality scores based on genotyping replicates;
- smart-thinning and final filtering;
- quality assessment based on replicates.

Also included are walkthroughs for analysis (2brad_denovo_analysis.txt):
- computing Weir and Cockerham Fst
- BayeScan
- ADMIXTURE
- fastSTRUCTURE
- dadi

IMPORTANT: Do not sequence 2bRAD libraries on a HiSeq 4000 lane alone! Invariant bases (adaptor, restriction site) will be read very poorly, resulting in read trimming issues. Mix 20% of PhiX libraries with your 2bRAD samples to avoid this problem, or share lane with someone else's non-2bRAD samples. (If you already did, no worries - there is a solution; contact me.)