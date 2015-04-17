Whole genome de novo genotyping with 2bRAD
------------------------------------------

Mikhail Matz, matz@utexas.edu

2bRAD has been described in Wang et al 2012 
http://www.nature.com/nmeth/journal/v9/n8/abs/nmeth.2023.html 

The advanced features of 2bRAD include removal of PCR duplicates, as well as advanced filtering and assessment of overall genotyping quality based on genotyping replicates

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

