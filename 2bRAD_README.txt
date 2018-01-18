2bRAD de novo and reference-based walkthrough
November 2, 2017
Mikhail Matz (matz@utexas.edu) - ask me if anything does not work

=============================================

INSTALLATIONS (you can skip these until needed):

------ vcftools: 

git clone https://github.com/vcftools/vcftools.git 
./autogen.sh
./configure --prefix=$HOME/bin/vcftools
make
make install

# adding the program directory to $PATH
cd
nano .bashrc
	export PATH=$HOME/bin/vcftools/bin:$PATH
	# press ctl-O, Enter, ctl-X

re-login to make PATH changes take effect

------- Moments: 

cd
git clone https://bitbucket.org/simongravel/moments.git 
cd moments
python setup.py build_ext --inplace

# add this to .bashrc, section 2:
  export PYTHONPATH=$PYTHONPATH:$HOME/moments
# re-login

cds
cd RAD

------- ANGSD: 

# install xz first from https://tukaani.org/xz/
cd
wget https://tukaani.org/xz/xz-5.2.3.tar.gz --no-check-certificate
tar vxf xz-5.2.3.tar.gz 
cd xz-5.2.3/
./configure --prefix=$HOME/xz-5.2.3/
make
make install

# edit .bashrc:
nano .bashrc
   export LD_LIBRARY_PATH=$HOME/xz-5.2.3/lib:$LD_LIBRARY_PATH
   export LIBRARY_PATH=$HOME/xz-5.2.3/lib:$LIBRARY_PATH
   export C_INCLUDE_PATH=$HOME/xz-5.2.3/include:$C_INCLUDE_PATH
logout
# re-login

# now, install htslib:
cd
git clone https://github.com/samtools/htslib.git
cd htslib
make CFLAGS=" -g -Wall -O2 -D_GNU_SOURCE -I$HOME/xz-5.2.3/include"

cd
git clone https://github.com/ANGSD/angsd.git 
cd angsd
make HTSSRC=../htslib

# now adding ANGSD to $PATH
cd
nano .bashrc
# section 2:
   export PATH=$HOME/angsd:$PATH
   export PATH=$HOME/angsd/misc:$PATH
# save (Ctl-O, Ctl-X)

-------  ngsTools (incl. ngsCovar) :

cd ~/bin
git clone https://github.com/mfumagalli/ngsPopGen.git
cd ngsPopGen
make
mv ngs* ..
cd -

-------  NGSadmix :
cd ~/bin/
wget popgen.dk/software/download/NGSadmix/ngsadmix32.cpp 
g++ ngsadmix32.cpp -O3 -lpthread -lz -o NGSadmix
cd -

-------  stairwayPlot :

# project page: https://sites.google.com/site/jpopgen/stairway-plot
cdw
# get version from June 2016 (v2beta)
wget https://www.dropbox.com/s/tj4i02n36abwjl6/stairway_plot_v0.2.zip?dl=0
mv file stairway_plot_v0.2.zip
unzip stairway_plot_v0.2.zip

-------  ADMIXTURE
cd ~/bin/
wget https://www.genetics.ucla.edu/software/admixture/binaries/admixture_linux-1.23.tar.gz --no-check-certificate
tar vxf admixture_linux-1.23.tar.gz 
mv admixture_linux-1.23/admixture .
cd -

-------  plink 1.9:
cd ~/bin
wget https://www.cog-genomics.org/static/bin/plink171103/plink_linux_x86_64.zip
unzip plink_linux_x86_64.zip
cd -

==============================================

# downloading and installing all 2bRAD scripts in $HOME/bin (or change to whatever directory you want)
cd
mkdir bin 
cd ~/bin 
# cloning github repositories
git clone https://github.com/z0on/2bRAD_denovo.git
git clone https://github.com/z0on/2bRAD_GATK.git
# move scripts to ~/bin from sub-directories
mv 2bRAD_denovo/* . 
mv 2bRAD_GATK/* . 
# remove now-empty directory
rm -rf 2bRAD_denovo 
rm -rf 2bRAD_GATK 

# designating all .pl and .py files (perl and python scripts) as executable
chmod +x *.pl 
chmod +x *.py
chmod +x *.R

# adding ~/bin to your $PATH 
cd
nano .bashrc
# paste this where appropriate (note: .bashrc configuration might be specific to your cluster, consult your sysadmin if in doubt)
   export PATH=$HOME/bin:$PATH

# Ctl-o, Ctl-x  (to save and exit in nano)
# log out and re-login to make sure .bashrc changes took effect

# does it work?
# try running a script from $HOME:
cd
2bRAD_trim_launch.pl
# if you get "command not found" something is wrong 

# ==============================================
#      Genome reference placement (if you have it)

# will need bowtie2, samtools, and picard. They are pre-installed as modules on TACC; you will have to install them if you don't have these modules on your cluster. 
module load perl
module load bowtie
module load samtools
module load picard-tools

# assuming we have a fasta file mygenome.fasta and it lives in the directory $WORK/db

export GENOME_FASTA=$WORK/db/mygenome.fasta
export GENOME_DICT=$WORK/db/mygenome.dict 

# indexing genome for bowtie2 mapper
bowtie2-build $GENOME_FASTA $GENOME_FASTA

samtools faidx $GENOME_FASTA

export GENOME_DICT=$WORK/db/mygenome.dict 
java -jar $TACC_PICARD_DIR/picard.jar CreateSequenceDictionary R=$GENOME_FASTA  O=$GENOME_DICT

#==================
# Step 1: Splitting by in-read barcode, deduplicating and quality-filtering the reads

# creating a file of commands to run (assuming reads are in fastq files, one file per sample.
# (if samples were spread across several lanes, concatenate them first using 
ngs_concat.pl) 
2bRAD_trim_launch_dedup.pl fastq > trims

# Note: use this command instead of the one above if you have 2bRAD libraries without degenerate tag but with 4-base in-line barcode:
2bRAD_trim_launch.pl fastq barcode2=4 > trims
# And if you have the original-style libraries without degenerate tag and without in-line barcode, use this:
2bRAD_trim_launch.pl fastq > trims
# AND IF you sequenced your double-barcoded, deduplicatable libraries on HiSeq 4000 alone
# (resulting in poor quality at restriction site and adaptor bases) and you have used BcgI enzyme, use this:
2bRAD_trim_launch_dedup2.pl fastq > trims

# the file trims now contains a long list of commands, one per fastq file, that we will need to execute. TACC uses Launcher module to enable parallele launch of multiple commands: https://github.com/TACC/launcher 
# I very highly recommend it - it should work on any SLURM-based system. Otherwise, if you don't know how to run hundreds of parallel jobs on your cluster, consult your IT support - there must be a way. 
# If push comes to shove, you can always execute all these commands serially (one by one), by saying 
bash trims

# do we have expected number of *.tr0 files created?
ls -l *.tr0 | wc -l

# quality filtering using fastx_toolkit (install fastx_toolkit if you don't have this module)
module load fastx_toolkit
ls *.tr0 | perl -pe 's/^(\S+)\.tr0$/cat $1\.tr0 \| fastq_quality_filter -q 20 -p 100 >$1\.trim/' >filt0

# NOTE: run the next line ONLY if your qualities are 33-based 
# (if you don't know just try to see if it works eventually, if you get errors from fastx_toolkit, try the other one):
	cat filt0 | perl -pe 's/filter /filter -Q33 /' > filt
#if you did NOT run the line above, run this one (after removing # symbol):
#	mv filt0 filt

# execute all commands in filt file (serial or parallel using Launcher, if your system allows) 

# do we have expected number of *.trim files created?
ls -l *.trim | wc -l

#====================
# denovo RAD business (skip to next "#========" if you have genome reference):

# 'uniquing' ('stacking') individual trimmed fastq reads:
ls *.trim | perl -pe 's/^(.+)$/uniquerOne.pl $1 >$1\.uni/' >unii

# execute all commands written to unii...

# Done! do you have .uni for all your samples?... 
ls -l *.uni | wc -l  

# merging uniqued files, creating mergedUniqTags.fasta for clustering 
# (set minDP and especially minInd higher if you want to crank up filtering against rare and possibly erroneous tag sequences)
mergeUniq.pl uni minDP=2 minInd=2 >mydataMerged.uniq

# Done! how many reads went into analysis?
tail mer.e*

# discarding tags that have more than 7 observations without reverse-complement
head -1 all.uniq >all.tab
awk '!($3>7 && $4==0)' all.uniq >>all.tab

# creating fasta file out of merged and filtered tags:
awk '{print ">"$1"\n"$2}' all.tab | tail -n +3 > all.fasta

# clustering reads into loci using cd-hit
# clustering allowing for up to 3 mismatches (-c 0.91); the most abundant sequence becomes reference
cd-hit-est -i all.fasta -o cdh_alltags.fas -aL 1 -aS 1 -g 1 -c 0.91 -M 0 -T 0  

#------------
# making fake reference genome (of 30 chromosomes) out ot cd-hit cluster representatives
# need bowtie2, samtools and picard_tools for indexing

concatFasta.pl fasta=cdh_alltags.fas num=30

# formatting fake genome
export GENOME_FASTA=cdh_alltags_cc.fasta
export GENOME_DICT=cdh_alltags_cc.dict 
module load perl
module load bowtie
bowtie2-build $GENOME_FASTA $GENOME_FASTA
module load samtools
samtools faidx $GENOME_FASTA
module load picard-tools
java -jar $TACC_PICARD_DIR/picard.jar CreateSequenceDictionary R=$GENOME_FASTA  O=$GENOME_DICT

#==============
# Mapping reads to reference (reads-derived fake one, or real) and formatting bam files 

# for denovo: map with bowtie2 with end-to-end matching 
export GENOME_FASTA=cdh_alltags_cc.fasta
2bRAD_bowtie2_launch.pl '\.trim$' $GENOME_FASTA | perl -pe 's/--score-min L,16,1 --local -L 16 //'> bt2

# for reference-based: mapping with --local option, enables clipping of mismatching ends (guards against deletions near ends of RAD tags)
export GENOME_FASTA=mygenome.fasta
2bRAD_bowtie2_launch.pl '\.trim$' $GENOME_FASTA > bt2

# execute all commands written to bt2...

# what are mapping efficiencies? 
cat maps.e*

# find out mapping efficiency for a particular input file (O9.fastq in this case)
# (assuming all input files have different numbers of reads)
grep -E '^[ATGCN]+$' O9.*trim | wc -l | grep -f - maps.e* -A 4 

ls *.bt2.sam > sams
cat sams | wc -l  # number should match number of trim files

# next stage is compressing, sorting and indexing the SAM files, so they become BAM files:
cat sams | perl -pe 's/(\S+)\.sam/samtools import \$GENOME_FASTA $1\.sam $1\.unsorted\.bam && samtools sort -o $1\.sorted\.bam $1\.unsorted\.bam && java -Xmx5g -jar \$TACC_PICARD_DIR\/picard\.jar AddOrReplaceReadGroups INPUT=$1\.sorted\.bam OUTPUT=$1\.bam RGID=group1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$1 && samtools index $1\.bam/' >s2b

rm *sorted*
ls *bam | wc -l  # should be the same number as number of trim files

# BAM files are the input into various genotype calling / popgen programs, this is the main interim result of the analysis. Archive them.

#===================== A  N  G  S  D =====================

# "FUZZY genotyping" with ANGSD - without calling actual genotypes but working with genotype likelihoods at each SNP. Optimal for low-coverage data (<10x).
# if your coverage is >10x, go to GATK section below

# install ANGSD first (see Installations section above)

# listing all bam filenames 
ls *bam >bams

#----------- assessing base qualities and coverage depth

export GENOME_REF=mygenome.fasta

# angsd settings:
# -minMapQ 20 : only highly unique mappings (prob of erroneous mapping = 1%)
# -baq 1 : realign around indels (not terribly relevant for 2bRAD reads mapped with --local option) 
# -maxDepth : highest total depth (sum over all samples) to assess; set to 10x number of samples

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -baq 1 -ref $GENOME_REF -maxDepth 1000"

# T O   D O : 
TODO="-doQsDist 1 -doDepth 1 -doCounts 1"

# in the following line, -r argument is one chromosome or contig to work with (no need to do this for whole genome as long as the chosen chromosome or contig is long enough)
# (look up lengths of your contigs in the header of *.sam files)
angsd -b bams -r chr1 -GL 1 $FILTERS $TODO -P 1 -out dd 

# summarizing results (using modified script by Matteo Fumagalli)
Rscript ~/bin/plotQC.R dd  
cat dd.info 
# scp dd.pdf to laptop to look at distribution of base quality scores, fraction of sites in each sample passing coverage thresholds, and fraction of sites passing genotyping rates cutoffs. Use these to guide choices of -minQ,  -minIndDepth and -minInd filters in subsequent ANGSD runs

#--------------- population structure

# Note: PCA and Admixture are not supposed to be run on data that contain clones (or genotyping replicates); manually remove them from bams list. If you want to detect clones, however, do keep the replicates and analyse identity-by-state (IBS) matrix (explained below)

# Generating genotype likelihoods from highly confident (non-sequencing-error) SNPs

# F I L T E R S :
# (ALWAYS record your filter settings and explore different combinations to confirm that results are robust. )
# Suggested filters :
# -minMapQ 20 : only highly unique mappings
# -minQ 30 : only highly confident base calls
# -minInd 50 : the site must be genotyped in at least 50 individuals (note: set this to at least 80% of your total number of your individuals)
# -snp_pval 1e-5 : high confidence that the SNP is not just sequencing error 
# -minMaf 0.05 : only common SNPs, with allele frequency 0.05 or more. Consider raising this to 0.1 for population structure analysis.
# Note: the last two filters are very efficient against sequencing errors but introduce bias against true rare alleles. It is OK (and even desirable) - UNLESS we want to do AFS analysis. We will generate data for AFS analysis in the next part.
# also adding filters against very badly non-HWE sites (such as, all calls are heterozygotes => lumped paralog situation) and sites with really bad strand bias:
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 30 -minInd 50 -snp_pval 1e-5 -minMaf 0.05 -dosnpstat 1 -doHWE 1 -hwe_pval 1e-5 -sb_pval 1e-5"

# THINGS FOR ANGSD TO DO : 
# -GL 1 : samtools likelihood model
# -doGlf 2 : output beagle format (for admixture)
# -doPost 1 : output posterior allele frequencies based on HWE prior
# -doGeno 32 : binary genotype likelihoods format (for ngsCovar => PCA)
# -doMajorMinor 1 : infer major and minor alleles from data (not from reference)
# -makeMatrix 1 -doIBS 1 -doCov 1 : identity-by-state and covariance matrices based on single-read resampling (robust to variation in coverage across samples)
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 32 -doPost 1 -doGlf 2"

# Starting angsd with -P the number of parallel processes. Funny but in many cases angsd runs faster on -P 1
angsd -b bams -GL 1 $FILTERS $TODO -P 1 -out myresult

# how many SNPs?
NSITES=`zcat myresult.beagle.gz | wc -l`
echo $NSITES

# covariance matrix for PCA (if coverage is approximately equal across samples; set nind to your number of samples):
gunzip myresult.geno.gz
ngsCovar -probfile myresult.geno -outfile myresult.covar -nind 61 -nsites $NSITES -call 0 -norm 0 

# if coverage is highly unequal among samples, use myresult.covMat and myresult.ibsMat from angsd run for PCoA and PCA. In fact, using these results would be conservative in any case.

# ADMIXTURE for K from 2 to 5
for K in `seq 2 5` ; 
do 
NGSadmix -likes myresult.beagle.gz -K $K -P 10 -o mydata_k${K};
done

# scp *Mat, *covar, *qopt and bams files to laptop, use angsd_ibs_pca.R to plot PCA and admixturePlotting_v4.R to plot ADMIXTURE

#==========================
# ANDSD => SFS for demographic analysis

# make separate files listing bams for each population (without clones and replicates)
# assume we have two populations, pop0 and pop1, 20 individuals each, with corresponding bams listed in pop0.bams and pop1.bams
# estimating site frequency likelihoods for each population - no filters except by genotyping rate (minInd, minIndDepth)
export GENOME_REF=mygenome.fasta
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 30 -minIndDepth 2"
TODO="-doSaf 1 -anc $GENOME_REF"
# In the following lines, set minInd to 80-90% of 2x pop sample size; i.e., if sample size is 20 set to 2*20*0.9: 36)
angsd -b pop0.bams -GL 1 -P 1 -minInd 36 $FILTERS $TODO -out pop0
angsd -b pop1.bams -GL 1 -P 1 -minInd 36 $FILTERS $TODO -out pop1

# generating per-population SFS
realSFS pop0.saf.idx >pop0.sfs
realSFS pop1.saf.idx >pop1.sfs

# generating dadi-like posterior counts based on sfs priors
realSFS dadi pop0.saf.idx pop1.saf.idx -sfs pop0.sfs -sfs pop1.sfs -ref $GENOME_REF -anc $GENOME_REF >dadiout

# converting to dadi-snp format understood by dadi an Moments:
# (numbers after the input file name are numbers of individuals sampled per population)
realsfs2dadi.pl dadiout 20 20 >2pops_dadi.data

#=====================
# 2d AFS analysis using Moments

# install Moments (python package, see the beginning of this file)

# get Misha's moments scripts collection
git clone https://github.com/z0on/AFS-analysis-with-moments.git

# print folded 2d SFS - for denovo or when mapping to genome of the studied species
# (change numbers to 2 x 0.9 x number of samples for in each pop):
2dAFS_fold.py 2pops_dadi.data pop0 pop1 36 26

# print unfolded 2d SFS - if mapping to genome of sister species
# (change numbers to 2 x 0.9 x number of samples for in each pop):
2dAFS.py 2pops_dadi.data pop0 pop1 36 36

# S2M model for pop0 and pop1: ancestral population splits into two different constant sizes, with asymmetrical migration 
# numbers are 2x number of samples for in each pop, and then starting values of parameters: nu1, nu2, T, m12, m21, and fraction of misidentified ancestral states (for unfolded). Parameters will be perturbed each time so all runs will be different
# run this same command 20 or more times to ensure convergence: 
S2M.py 2pops_dadi.data pop0 pop1 36 36 1 1 1 1 1 0.01
# folded (for denovo):
S2M_fold.py 2pops_dadi.data pop0 pop1 36 36 1 1 1 1 1

# IM2 model: ancestral population splits into two different sizes which then experience different exponential growth/decline, with asymmetrical migration 
# Folded IM2 model for pop0 and pop1 (parameters: nu1_0,nu2_0,nu1,nu2,T,m12,m21):
IM2_fold.py 2pops_dadi.data pop0 pop1 36 36 1 1 1 1 1 1 1

#==========================

#        G  A  T  K

# ("hard-call" genotyping, use only for high-coverage data, >10x after deduplication)

# invoke modules (skip if using your own installations)
module load gatk
module load picard-tools

# UNLESS you are working on TACC, edit these accordingly and execute:
export TACC_GATK_DIR=/where/gatk/is/installed/
export TACC_PICARD_DIR=/where/picard/is/installed/
# NOTE that you will have to execute the above two lines every time you re-login (put them into your .bashrc)

# (see above about formatting of the reference genome - same thing here although it is a real, not reads-derived genome)
export GENOME_REF=mygenome.fasta
ls *.bam > bams

# writing command script
echo '#!/bin/bash
#SBATCH -J gt
#SBATCH -n 40
#SBATCH -N 1
#SBATCH -p development
#SBATCH -o gt.o%j
#SBATCH -e gt.e%j
#SBATCH -t 2:00:00
#SBATCH -A mega2014
#SBATCH --mail-type=ALL
#SBATCH --mail-user=matz@utexas.edu
java -jar $TACC_GATK_DIR/GenomeAnalysisTK.jar -T UnifiedGenotyper \
-R $GENOME_REF -nt 40 -nct 1 \
--genotype_likelihoods_model SNP \' >unig2
cat bams | perl -pe 's/(\S+\.bam)/-I $1 \\/' >> unig2
echo '-o primary.vcf ' >> unig2
bash unig2

#----------
# Variant quality score recalibration (VQSR)

# making a tab-delimited table of clone (replicate) sample pairs
nano clonepairs.tab
# paste names of bams that are clone pairs, tab delimited, one pair per line; for example
# Ctl-O , enter, Ctl-X

# extracting "true snps" subset (reproducible across replicates) 
# parameter hetPairs can vary depending on replication scheme (3 is good when you have triplicates)
replicatesMatch.pl vcf=primary.vcf replicates=clonepairs.tab hetPairs=2 max.het=0.5 > vqsr.vcf

# determining transition-transversion ratio for true snps (will need it for tranche calibration)
vcftools --vcf vqsr.vcf --TsTv-summary
# Ts/Tv ratio: 1.44  # put your actual number into the next code chunk, --target_titv

# creating recalibration models
export GENOME_REF=mygenome.fasta
java -jar $TACC_GATK_DIR/GenomeAnalysisTK.jar -T VariantRecalibrator \
-R $GENOME_REF -input primary.vcf -nt 12 \
-resource:repmatch,known=true,training=true,truth=true,prior=30  vqsr.vcf \
-an QD -an MQ -an FS -mode SNP --maxGaussians 6 \
--target_titv 1.44 -tranche 85.0 -tranche 90.0 -tranche 95.0 -tranche 99.0 -tranche 100 \
-recalFile primary.recal -tranchesFile recalibrate.tranches -rscriptFile recalibrate_plots.R 

# examine output and recalibrate*.pdf files - see which tranche to choose the one before TsTv dropoff
# next chunk assumes we are choosing tranche 95

# applying recalibration (95% tranche)
export GENOME_REF=mygenome.fasta
java -jar $TACC_GATK_DIR/GenomeAnalysisTK.jar -T ApplyRecalibration \
-R $GENOME_REF -input primary_n.vcf -nt 12 \
--ts_filter_level 95.0 -mode SNP \
-recalFile primary.recal -tranchesFile recalibrate.tranches -o recal.vcf

#---------------
# Applying filters

# identifying poorly genotyped individuals
vcftools --vcf recal.vcf --het
# look at number of sites genotyped per individual (4th column): 
cat out.het 
# see if some samples are much lower in the number of sites than others
# for example, if you want to remove samples showing less than 40000 sites:
cat out.het | awk '$4<40000' | cut -f 1  > underSequenced
cat underSequenced

# applying filter and selecting polymorphic biallelic loci genotyped in 90% or more individuals
# (harsh genotyping rate cutoff is strongly recommended for best quality and to avoid RAD loci affected by null 
# alleles because of mutations in restriction site)
vcftools --vcf recal.vcf --remove underSequenced --remove-filtered-all --max-missing 0.9  --min-alleles 2 --max-alleles 2 --recode-INFO-all --recode --out filt

# selecting only polymorphic sites (they all are in denovo pipeline!) and sites with no excess heterozygosity
grep -E "#|0/1|0/0.+1/1|1/1.+0/0" filt.recode.vcf >polymorphs.vcf
hetfilter.pl vcf=polymorphs.vcf maxhet=0.5 >best.vcf

#---------------
# Final touches

# genotypic match between pairs of replicates (the most telling one is the last one, HetsDiscoveryRate - fraction of correctly called heterozygotes; if it is under 90% perhaps use fuzzy genotyping with ANGSD - see above)	
repMatchStats.pl vcf=best.vcf replicates=clonepairs.tab 

# looking at per-individual inbreeding 
# positive - excess false homozygotes (due to poor coverage); negative - false heterozygotes (possibly lumped paralogs)
vcftools --vcf best.vcf --het
cat out.het

# create a file listing clones (and low-site/high-homozygosity individuals, if any) to remove
cat clonepairs.tab | cut -f 2 >clones2remove

# removing clones and badly genotyped ones
vcftools --vcf best.vcf --remove clones2remove --recode --recode-INFO-all --out final

# thinning for Fst / PCA / ADMIXTURE  (choosing one SNP per tag with max allele frequency):
thinner.pl vcf=final.recode.vcf criterion=maxAF >thinMaxaf.vcf

#----------
# ADMIXTURE using hard-call data (vcf)

# for denovo: creating a dataset with fake "chr" chromosome designations
cat thinMaxaf.vcf | perl -pe 's/tag(\d)(\d+)\t(\d+)/chr$1\t$3$2/'>thinMaxChrom.vcf

# reformatting VCF into plink binary BED format (install plink first! see beginning of this readme)
plink --vcf thinMaxChrom.vcf --make-bed --out rads

# ADMIXTURE with cross-validation to select K  (install it first, see beginning of this readme)
for K in 1 2 3 4 5; do admixture --cv rads.bed $K | tee log${K}.out; done

# minimal cross-validation error = optimal K
grep -h CV log*.out

# create a two-column tab-delimited file inds2pops: bams listed in column 1, their pop designations - in column 2
# scp the *.Q and inds2pops to laptop, plot it in R:
tbl=read.table("rads.3.Q")
barplot(t(as.matrix(tbl)), col=rainbow(5),xlab="Individual #", ylab="Ancestry", border=NA)

# or, more fancy, use admixturePlotting_V4.R

#---------------
# SFS analysis of hard-call (vcf) data

# create a two-column tab-delimited file inds2pops: bams listed in column 1, their pop designations - in column 2

# converting vcf to dadi-snp format (the output is file with _dadi.data extension)
vcf2dadi.pl final.recode.vcf inds2pops

# use Moments scripts (see above, just before GATK section)









