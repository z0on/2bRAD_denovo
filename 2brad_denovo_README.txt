2bRAD de novo walkthrough
September 15, 2014 Mikhail Matz (matz@utexas.edu)

The idea is to copy the chunks separated by empty lines below and paste them into your cluster 
terminal window consecutively. 

The lines beginning with hash marks (#) are explanations and additional instructions - 
please make sure to read them before copy-pasting. 

In addition to the scripts coming with this distribution,
you will need the following software installed and available 
(note: TACC already has them as modules):

python: http://www.python.org/getit/
fastx_toolkit: http://hannonlab.cshl.edu/fastx_toolkit/download.html
cd-hit: https://cdhit.googlecode.com/files/cd-hit-v4.6.1-2012-08-27.tgz
vcftools: http://vcftools.sourceforge.net/ 

Occasionally, an example screen output is listed after the program calls below.

NOTE: this walkthrough has been written for lonestar cluster of the Texas
Advanced Computer Center, which has 12 cores per node and uses Sun Grid Engine 
(SGE) batch processing system. To adopt the walkthough to your cluster you 
would need to edit the launcher_creator.py script to make its default settings
compatible with your cluster. 

==============================================

# unzipping multiple *.gz files
ls *.gz | perl -pe 's/^(.+)$/gunzip $1/' >gunz
launcher_creator.py -j gunz -n gunz -l gj
qsub gj

# (ecogeno2013: skip this one)
# concatenating the reads files (make sure to unzip them first!):
# NOTE: Only run this chunk if your samples are spread over several lanes.
# Assuming your filed have the extension fastq (edit that in the line below if not),
# replace "2b_(.+)_L00" below with the actual pattern to recognize sample
# identifier in the file name. The pattern above will work for file names such as
# Sample_2b_K208_L007_R1.cat.fastq (the identifier is K208)
ngs_concat.pl fastq "2b_(.+)_L00"

# Trimming the reads, splitting them by secondary barcode and removing PCR duplicates
# The sampleID parameter here, '3', defines a chunk in the filename (separated
# by underscores) that is to be kept; 
# for example, for a filename like this: Sample_2b_M11_L006_R1.cat.fastq 
# it would be reasonable to specify 'sampleID=3' to keep only 'M11' as a 
# sample identifier. If the sampleID parameter is not specified, the whole filename up to the first dot will be kept.
# NOTE: if you ran ngs_concat.pl (above), run not the next but the second-next line (after removing # symbol) :
2bRAD_trim_launch_dedup.pl fastq > trims
launcher_creator.py -j trims -n trims -l trimjob
qsub trimjob

# NB: use this command instead of the one above if you have 2bRAD libraries without 
# degenerate tag but with 4-base in-line barcode:
2bRAD_trim_launch.pl fastq barcode2=4 > trims
# And if you have the original-style libraries without degenerate tag and without in-line barcode, use this:
2bRAD_trim_launch.pl fastq > trims


# quality filtering using fastx_toolkit
module load fastx_toolkit
ls *.tr0 | perl -pe 's/^(\S+)\.tr0$/cat $1\.tr0 \| fastq_quality_filter -q 20 -p 90 >$1\.trim/' >filt0

# NOTE: run the next line ONLY if your qualities are 33-based 
# (if you don't know just try to see if it works eventually, if you get errors from fastx_toolkit, try the other one):
	cat filt0 | perl -pe 's/filter /filter -Q33 /' > filt
#if you did NOT run the line above, run this one (after removing # symbol):
#	mv filt0 filt

launcher_creator.py -j filt -n filt -l filtjob
qsub filtjob


# uniquing individual trimmed fastq files:
ls *.trim | perl -pe 's/^(.+)$/uniquerOne.pl $1 >$1\.uni/' > unii
launcher_creator.py -j unii -n uni1 -l uni1
qsub uni1

ll *.uni | wc -l  
# do you have .uni for all your samples?... If not, rerun the chunk above

# merging uniqued files, creating mergedUniqTags.fasta for clustering :
echo "mergeUniq.pl uni >mydataMerged.uniq" > mer 
launcher_creator.py -j mer -n mjob -l mjob
cat mjob | perl -pe 's/12way 12/1way 12/' >mjobb
qsub mjobb

tail mjob.e*
# considering reads seen at least 5 times:
# 1014968293 reads processed

# clustering allowing for up to 3 mismatches (-c 0.91); the most abundant sequence becomes reference
echo 'cd-hit-est -i mergedUniqTags.fasta -o cdh_alltags.fas -aL 1 -aS 1 -g 1 -c 0.91 -M 0 -T 0' >cdh
launcher_creator.py -j cdh -n cdhjobB -l cdhjob
cat cdhjob | perl -pe 's/12way 12/1way 12/' >cdhjobB
qsub cdhjobB


# assembling locus(cluster)-annotated table, filtering by read depth (def:5) and strand bias (def:5)
echo 'uniq2loci.pl mydataMerged.uniq cdh_alltags.fas.clstr cdh_alltags.fas >cdh_alltags.ul' >uloc
launcher_creator.py -j uloc -n uloc -l ulocj
cat ulocj | perl -pe 's/12way 12/1way 12/' >ulocjob
qsub ulocjob

cat uloc.e*
# 2290109	total tags
# 2290109	pass depth cutoff 5
# 1580389	pass strand bias cutoff 5


# calling genotypes while applying mild total depth (aobs, default 10), 
# allele bias (abias, def 10) and strand bias (strbias, def 10) filters
# also requiring an allele to appear in at least 2 individuals (ind) 
# maximum allowed heterozygosity at the locus: def 0.5 (hetero)
# a homozygous genotype will not be called if read depth at a locus in a sample is less than 3 (mindp, default 3)
# (all these parameters are adjustable)
echo 'haplocall_denovo.pl cdh_alltags.ul' >hcdn
launcher_creator.py -j hcdn -n hcdn -l hcdnj
cat hcdnj | perl -pe 's/12way 12/1way 12/' >hcdnjob
qsub hcdnjob

tail -100 hcdn.o*
#	Allele filters:
#	1580060	raw alleles
#	1083572	with 10 or more reads
#	1019901	in 2 or more samples
#	971260	pass strand bias cutoff 10
#	698160	pass allele bias cutoff 10
#	--------------------
#	Locus filters:
#	1064910	total
#	449701	remain after applying allele filters
#	437460	have less than 0.5 fraction of heterozygotes
#	416586 	with not more than 2 alleles
#	92902	genotyped in 25% of samples
#	43641	polymorphic


# making a tab-delimited table of clone (replicate) sample pairs
nano clonepairs.tab
# ecogeno2013: paste this :
# everyone else: edit accordingly 
K210	K212
K212	K213
K213	K216
K211	K219

# IF YOU HAVE REPLICATES (strongly recommended):
# extracting "true snps" subset (reproducible across replicates)
replicatesMatch.pl vcf=cdh_alltags.ul_Variants_count10_ab10_sb10_clip0.vcf replicates=clonepairs.tab polyonly=1 > vqsr.denovo.vcf
#	76865 total SNPs
#	27874 passed
#	8642 alt alleles
#	8616 polymorphic
#	1067 written

# same for haplotypes:
replicatesMatch.pl vcf=cdh_alltags.ul_Vhap_count10_ab10_sb10_clip0.vcf replicates=clonepairs.tab polyonly=1 > vqsr.vhap.vcf

# non-parametric recalibration (strand bias not used, -nosb)
# the goal is to generate a composite filter and determine the setting that gives maximum "gain"	
recalibrateSNPs.pl vcf=cdh_alltags.ul_Variants_count10_ab10_sb10_clip0.vcf true=vqsr.denovo.vcf -nosb -notp >denovo.recal.vcf
#------------------------
# 8.01%	at qual <1 (7.01% gain)
# 25.86%	at qual <5 (20.86% gain)
# 31.10%	at qual <10 (21.10% gain)
# 35.53%	at qual <15 (20.53% gain)
# 42.87%	at qual <20 (22.87% gain)
# 49.93%	at qual <30 (19.93% gain)
#------------------------

# also recalibrating haplotype-wise calls
recalibrateSNPs.pl vcf=cdh_alltags.ul_Vhap_count10_ab10_sb10_clip0.vcf true=vqsr.vhap.vcf -nosb -notp >denovo.vhap.recal.vcf
# ------------------------
# 7.57%	at qual <1 (6.57% gain)
# 30.89%	at qual <5 (25.89% gain)
# 35.85%	at qual <10 (25.85% gain)
# 42.40%	at qual <15 (27.40% gain)
# 45.72%	at qual <20 (25.72% gain)
# 53.07%	at qual <30 (23.07% gain)
# ------------------------


# if you DO NOT have replicates:
# you will have to choose filtering criteria based solely on quantiles observed in the your
# dataset. To determine the distribution of AB, SB and DP:
filterStats.pl vcf=alltags.ul_Variants_count10_ab10_sb10_clip0.vcf
# AB quantiles:
# 99%:	>3
# 95%:	>15
# 90%:	>21
# 85%:	>26
# 80%:	>29
# 70%:	>35
# 60%:	>39
# 50%:	>43
# 40%:	>46
# 30%:	>49
# 20%:	>53
# 10%:	>60
# 
# SB quantiles:
# 99%:	>24
# 95%:	>34
# 90%:	>41
# 85%:	>47
# 80%:	>52
# 70%:	>61
# 60%:	>68
# 50%:	>74
# 40%:	>78
# 30%:	>83
# 20%:	>87
# 10%:	>92
# 
# meanDP quantiles:
# 99%:	7 - 424
# 95%:	10 - 143
# 90%:	16 - 122
# 85%:	21 - 110
# 80%:	23 - 102
# 70%:	25 - 90
# 60%:	28 - 81
# 50%:	31 - 73
# 40%:	34 - 67
# 30%:	37 - 62
# 20%:	40 - 57
# 10%:	44 - 52

# I recommend using AB and DP filters (SB and TP don't seem to be too useful at this point)

# if you did replicates-based recalibration and your best gain was observed at 20:
# also, minimum genotype quality we want is Q15, biallelic only, present in 80% of samples:	
vcftools --vcf alltags.ul_Variants_count10_ab10_sb10_clip0.vcf --minQ 20 --minGQ 15 --min-alleles 2 --max-alleles 2 --max-missing 0.75 --recode --recode-INFO-all --out denovo.filt0

# if you did NOT do recalibration and want (recommended) to keep loci within best 80% for each filter (AB and meanDP):
# consult the output from filterStats.pl (above) to see which filter values correspond to 80%
# in vcftools, AB would correspond to --minQ, meanDP range is set by --min-meanDP and --max-meanDP:
vcftools --vcf alltags.ul_Variants_count10_ab10_sb10_clip0.vcf --minQ 29 --min-meanDP 23 --max-meanDP 102 --recode --recode-INFO-all --out denovo.filt0
# After filtering, kept 132 out of 132 Individuals
# After filtering, kept 26023 out of a possible 76865 Sites

# thinning SNP dataset (leaving one snp per tag, with the highest minor allele frequency)
thinner.pl vcf=denovo.filt0.recode.vcf > thinDenov.vcf

# evaluating accuracy (the most important one is Heterozygote discovery rate, last column) 
# based on replicates
repMatchStats.pl vcf=thinDenov.vcf replicates=clonepairs.tab 
#	pair	gtyped	match	[ 00	01	11 ]	HetMatch	HomoHetMismatch	HetNoCall	HetsDiscoveryRate
#K210:K212	7328	7169(97.8%)	 [78%	17%	5% ]	1200	139	5	0.94	
#K212:K213	7369	7117(96.6%)	 [78%	17%	5% ]	1202	179	4	0.93	
#K213:K216	7313	7026(96.1%)	 [78%	17%	5% ]	1192	129	12	0.94	
#K211:K219	7312	7084(96.9%)	 [78%	18%	5% ]	1264	222	2	0.92	

# creating a file of replicates and other poorly sequenced samples (low sites number/high-homozygosity,
# use vcftools --het to look those up) to remove:
echo 'K210
K212
K216
K219
K3
O4
K4
O5'>clones2remove

# creating final, thinned dataset without clones
vcftools --vcf thinDenov.vcf --remove clones2remove --recode --recode-INFO-all --out denovo.filt

#-------------------
# Now the same recalibration and filtering for haplotype-wise calls 
# (you can skip this unless you plan coalescent analysis)
vcftools --vcf denovo.vhap.recal.vcf --minQ 20 --minGQ 15 --max-missing 0.8 --recode --recode-INFO-all --out denovo.vhap.filt

# assessing quality
repMatchStats.pl vcf=denovo.vhap.filt.recode.vcf replicates=clonepairs.tab 

# making final set with no replicates
vcftools --vcf denovo.vhap.recal.vcf --remove clones2remove --min-alleles 2 --minQ 20 --minGQ 15 --max-missing 0.8 --recode --recode-INFO-all --out denovo.vhap.filt
# After filtering, kept 9251 out of a possible 45373 Sites

# if you have a genome draft to check the percentage of truly unique loci:
# (but then again, if you have a genome draft better use reference-based GATK pipeline!)
vhap2fasta.pl denovo.vhap.filt.recode.vcf > dnloci.fasta
bowtie2 -L 16 -x /work/01211/cmonstr/amil_genome/amil_genome_fold_c.fasta -U dnloci.fasta -f >dnloci.sam
#9251 reads; of these:
#  9251 (100.00%) were unpaired; of these:
#    620 (6.70%) aligned 0 times
#    7922 (85.63%) aligned exactly 1 time
#    709 (7.66%) aligned >1 times
#93.30% overall alignment rate

# optional...

# For coalescent-based analysis (MIGRATE-N, IMa) - creating haplotypes dataset while 
keeping monomorphic tags ('mono=keep') 
echo 'haplocall_denovo.pl alltags.ul mindp=10 mono=keep' >hcdnm
launcher_creator.py -j hcdnm -n hcdnm -l hcdnmj
cat hcdnmj | perl -pe 's/12way 12/1way 12/' >hcdnmjob
qsub hcdnmjob

#=================================================

# see 2brad_denovo_analysis.txt for some analysis options
