Whole genome de novo genotyping with 2bRAD
------------------------------------------

Mikhail Matz, matz@utexas.edu

Current sample preparation protocol (feel free to comment!): https://docs.google.com/document/d/1am7L_Pa5JQ4sSx0eT5j4vdNPy5FUAtMZRsJZ0Ar5g9U/edit?usp=sharing

2bRAD has been described in Wang et al 2012 
http://www.nature.com/nmeth/journal/v9/n8/abs/nmeth.2023.html 

2bRAD features a very simple library prep protocol with no intermediate purification stages. It involves triple barcode scheme, two being standard Illumina barcodes and one in-read ligated barcode for pooling libraries in 12-plexes midway through library prep to further minimize prep costs. 2bRAD generates 36-base tags that can be sequenced on either strand. In a full-representation version it results in one high-quality genotyped tag every 2.5-3.5 kb on average. With reduced-representation with modified ligated adapters, the number of genotyped tags can be reduced 4- or 16-fold, to save sequencing costs in applications not requiring dense genome sampling (such as population structure analysis or linkage mapping). 

Unique features of 2bRAD bioinformatics are removal of PCR duplicates based on degenerate tag (allowing for 128-fold dynamic range of sequencing depth per allele), use of direct-reverse strand ratio as a quality fltering parameter, and filtering and quality assessment based on genotyping replicates.

de novo 2bRAD pipeline is similar in ideology to STACKS, but takes advantage of the fact that 2bRAD sequences both strands. The assembled RAD loci are formatted as a fake genome, which is then used as a reference - from there on, the pipeline follows the reference-based steps. Since mid-2017 we favor ANGSD as our main data processing toolkit (http://www.popgen.dk/angsd/index.php/ANGSD). The 2bRAD_README.sh file suggests basic steps to analyze genotyping quality, investigate population subdivision and obtain unbiased SFS for demographic analysis (using dadi or Moments). We also keep instructions for  GATK UnifiedGenotyper and variant quality score recalibration; these would be preferred for high-coverage data (10x or better). 

This repository contains the lab protocol for sample preparation, as well as scripts and walkthrough (2bRAD_README.txt) for:
- read trimming, quality filtering, removal of PCR duplicates
- assembling RAD loci
- creating and formatting a fake genome out of RAD loci
- mapping original reads to fake genome
- exploring coverage depth and genotyping rates

For <10x coverage data:
- performing PCA and ADMIXTURE analysis based on "fuzzy" genotyping in ANGSD
- deriving site frequency spectra using ANGSD and analyzing them in Moments

For >10x coverage data:
- calling variants using GATK's UnifiedGenotyper
- variant quality score recalibration (VQSR) based on genotyping replicates
- final thinning and filtering
- quality assessment based on replicates
- deriving site frequency spectra

IMPORTANT: Do not sequence 2bRAD libraries on a HiSeq 4000 lane alone! Invariant bases (adaptor, restriction site) will be read very poorly, resulting in read trimming issues. Mix 20% of PhiX libraries with your 2bRAD samples to avoid this problem, or share lane with someone else's non-2bRAD samples. (If you already did, no worries - there is a solution; contact me.)
