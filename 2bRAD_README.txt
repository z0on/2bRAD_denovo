2bRAD de novo and reference-based walkthrough
June 2018
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

-------  PopGenTools (angsd wrapper and format converter) :

git clone https://github.com/CGRL-QB3-UCBerkeley/PopGenTools.git
ll PopGenTools/
cp PopGenTools/* ~/bin

----- PGDspider :

cd ~/bin
wget http://www.cmpg.unibe.ch/software/PGDSpider/PGDSpider_2.0.7.1.zip
unzip PGDSpider_2.0.7.1.zip
cd -

----- Bayescan :

cd ~/bin
wget http://cmpg.unibe.ch/software/BayeScan/files/BayeScan2.1.zip
unzip BayeScan2.1.zip
cp BayeScan2.1/binaries/BayeScan2.1_linux64bits bayescan
chmod +x bayescan
rm -r BayeScan*

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
java -jar $TACC_PICARD_T_DIR/picard.jar CreateSequenceDictionary R=$GENOME_FASTA  O=$GENOME_DICT

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
# set minInd to 25% of all your individuals 
# (assuming we want to assemble fake reference genome out of these tags, we only
# need to capture major alleles at each locus)
# set minDP to even higher value (total depth), something like 2 x minInd
mergeUniq.pl uni minDP=40 minInd=20 >all.uniq

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
java -jar $TACC_PICARD_T_DIR/picard.jar CreateSequenceDictionary R=$GENOME_FASTA  O=$GENOME_DICT

#==============
# Mapping reads to reference (reads-derived fake one, or real) and formatting bam files 

# for denovo: map reads to fake genome with bowtie2: 
GENOME_FASTA=cdh_alltags_cc.fasta
>maps
for F in `ls *.trim`; do
echo "bowtie2 --no-unal -x $GENOME_FASTA -U $F -S $F.sam">>maps
done

# for reference-based: mapping with --local option, enables clipping of mismatching ends (guards against deletions near ends of RAD tags)
GENOME_FASTA=mygenome.fasta
2bRAD_bowtie2_launch.pl '\.trim$' $GENOME_FASTA > maps

# execute all commands written to maps...

>alignmentRates
for F in `ls *trim`; do 
M=`grep -E '^[ATGCN]+$' $F | wc -l | grep -f - maps.e* -A 4 | tail -1 | perl -pe 's/maps\.e\d+-|% overall alignment rate//g'` ;
echo "$F.sam $M">>alignmentRates;
done

ls *.sam > sams
cat sams | wc -l  # number should match number of trim files

# next stage is compressing, sorting and indexing the SAM files, so they become BAM files:
export GENOME_REF=cdh_alltags_cc.fasta
cat sams | perl -pe 's/(\S+)\.sam/samtools import \$GENOME_REF $1\.sam $1\.unsorted\.bam && samtools sort -o $1\.sorted\.bam $1\.unsorted\.bam && java -Xmx5g -jar \$TACC_PICARD_T_DIR\/picard\.jar AddOrReplaceReadGroups INPUT=$1\.sorted\.bam OUTPUT=$1\.bam RGID=group1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$1 && samtools index $1\.bam/' >s2b

# run all commands listed in s2b file

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
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 30 -minInd 50 -snp_pval 1e-5 -minMaf 0.05 -dosnpstat 1 -doHWE 1 -sb_pval 1e-2 -hetbias_pval 1e-2 -skipTriallelic 1"

# THINGS FOR ANGSD TO DO : 
# -GL 1 : samtools likelihood model
# -doGlf 2 : output beagle format (for admixture)
# -doPost 1 : output posterior allele frequencies based on HWE prior
# -doGeno 32 : binary genotype likelihoods format (for ngsCovar => PCA)
# -doMajorMinor 1 : infer major and minor alleles from data (not from reference)
# -makeMatrix 1 -doIBS 1 -doCov 1 : identity-by-state and covariance matrices based on single-read resampling (robust to variation in coverage across samples)
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 32 -doVcf 1 -doPost 1 -doGlf 2"

# Starting angsd with -P the number of parallel processes. Funny but in many cases angsd runs faster on -P 1
angsd -b bams -GL 1 $FILTERS $TODO -P 1 -out myresult

# how many SNPs?
NSITES=`zcat myresult.beagle.gz | wc -l`
echo $NSITES

# covariance matrix (I don't see a point of this, IBS matrix is the same thing but better - not sensitive to variation in coverage):
gunzip myresult.geno.gz
ngsCovar -probfile myresult.geno -outfile myresult.covar -nind 61 -nsites $NSITES -call 0 -norm 0 

# if coverage is highly unequal among samples, use myresult.covMat and myresult.ibsMat from angsd run for PCoA and PCA. In fact, using these results would be conservative in any case.

# NgsAdmix for K from 2 to 5
for K in `seq 2 5` ; 
do 
NGSadmix -likes myresult.beagle.gz -K $K -P 10 -o mydata_k${K};
done

# alternatively, to use real ADMIXTURE on called SNPs:
gunzip myresult.vcf.gz
plink --vcf myresult.vcf --make-bed --out myresult
for K in `seq 1 5`; \
do admixture --cv myresult.bed $K | tee myresult_${K}.out; done

# which K is least CV error?
grep -h CV myresult_*.out

# scp the *.Q and inds2pops files to laptop, plot it in R:
# use admixturePlotting2a.R to plot (will require minor editing - population names)

# scp *Mat, *covar, *qopt and bams files to laptop, use angsd_ibs_pca.R to plot PCA and admixturePlotting_v4.R to plot ADMIXTURE

#==========================
# ANDSD => SFS for demographic analysis

# make separate files listing bams for each population (without clones and replicates)
# assume we have two populations, pop0 and pop1, 20 individuals each, with corresponding bams listed in pop0.bams and pop1.bams

# filtering sites to work on - use only filters that do not distort allele frequency
# set minInd to 75-90% of the total number fo individuals in the project
# if you are doing any other RAD than 2bRAD or GBS, remove '-sb_pval 1e-2' from FILTERS
cat pop0.bams pop1.bams > all.bams
FILTERS="-minMapQ 30 -minQ 35 -minInd 36 -doHWE 1 -sb_pval 1e-2 -hetbias_pval 1e-2 -skipTriallelic 1"
DOS="-doMajorMinor 1 -doMaf 1 -dosnpstat 1 -dogeno 3 -doPost 2"
angsd -b all.bams -GL 1 $FILTERS $DOS -P 1 -out sfsSites

# extracting and indexing list of sites to make SFS from 
# filtering out sites where heterozygote counts comprise more than 50% of all counts (likely lumped paralogs)
zcat sfsSites.snpStat.gz | awk '($3+$4+$5+$6)>0' | awk '($12+$13+$14+$15)/($3+$4+$5+$6)<0.5' | cut -f 1,2  >sites2do 
angsd sites index sites2do

# estimating site frequency likelihoods for each population 
export GENOME_REF=mygenome.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF""
# In the following lines, set minInd to 75-90% of each pop's sample size
angsd -sites sites2do -b pop0.bams -GL 1 -P 1 $TODO -minInd 18 -out pop0
angsd -sites sites2do -b pop1.bams -GL 1 -P 1 $TODO -minInd 18 -out pop1

# generating per-population SFS
realSFS pop0.saf.idx >pop0.sfs
realSFS pop1.saf.idx >pop1.sfs

# generating dadi-like posterior counts based on sfs priors
realSFS dadi pop0.saf.idx pop1.saf.idx -sfs pop0.sfs -sfs pop1.sfs -ref $GENOME_REF -anc $GENOME_REF >dadiout

# converting to dadi-snp format understood by dadi an Moments:
# (numbers after the input file name are numbers of individuals sampled per population)
realsfs2dadi.pl dadiout 20 20 >2pops_dadi.data

#===================== 2d AFS analysis using Moments

# get Misha's Moments scripts collection
git clone https://github.com/z0on/AFS-analysis-with-moments.git
# set your $PATH to include the AFS-analysis-with-moments directory

# print folded 2d SFS - for denovo or when mapping to genome of the studied species
# (change numbers to 2 x 0.9 x number of samples for in each pop):
2dAFS_fold.py 2pops_dadi.data pop0 pop1 36 36

# print unfolded 2d SFS - if mapping to genome of sister species
# (change numbers to 2 x 0.9 x number of samples for in each pop):
2dAFS.py 2pops_dadi.data pop0 pop1 36 36

# ------ multimodel inference: fit a diversity of 2-population models, then select the best one based on AIC.
# there are models with a period of exponential growth ("IM" models), models with one, two or three different size and/or migration rate epochs ("SC" models, including models with no migration in some epochs), models with symmetrical migration ("sm" models, in other cases migration is asymmetrical), models with two types of genomic loci ("genomic islands") introgressing at different rates ("i" models), and some fun combinations thereof.
# differences between models are summarized in excel table moments_multimodels.xls

# read about multimodel inference here: 
# https://pdfs.semanticscholar.org/a696/9a3b5720162eaa75deec3a607a9746dae95e.pdf

# this HAS to be parallelized - we need to fit ~100 models 10 times to make sure each model converges at its best fit at least once.

# input line: the last four numbers are:
# - projections (2 x 0.9 x number of samples) for in each pop;
# - mutation rate per gamete per generation
# - generation time, in thousand years
IN="2pops_dadi.data pop0 pop1 36 36 0.02 0.005"

# after setting the $IN variable, execute all commands listed in the text file "allmodels_unfolded" (if your alleles are polarized into ancestral and derived, for example by mapping to a sister species genome) or "allmodels_folded" 

# collecting results while fixing broken lines
cat *.mom | perl -pe 's/RESULT(.+)(\d)\n/RESULT$1$2/' |perl -pe 's/RESULT(.+)(\d)\n/RESULT$1$2/' |perl -pe 's/RESULT(.+)(\d)\n/RESULT$1$2/' | perl -pe 's/RESULT(.+)([\d\s])\n/RESULT$1$2/' | grep RESULT > mmods.res

# extracting likelihoods and parameter numbers for AIC:
cut -f 2,3,4,5 -d " " mmods.res >likes

# use R script deltaAIC_multimodels.R to find best-fitting model.
# using model ID number that the R script will list: 
# - examine the model's graphic output (*.pdf of actual and modeled SFS, and *.png of the model graph)
# - grep fitted model parameters and their SDs from *.mom files

# the order of parameters are listed in files unfolded_params and folded_params. Typically pop size parameters are first, then times, then migration rates, then the fraction of genomic "islands" (in "i"  models), then percentage of misidentified ancestral states (in unfolded models).

#=========== Fst: global, per site, and per gene

# writing down 2d-SFS priors
realSFS pop0.saf.idx pop1.saf.idx -P 24 > p01.sfs ; realSFS fst index pop0.saf.idx pop1.saf.idx -sfs p01.sfs -fstout p01 

# global Fst between populations
realSFS fst stats p01.fst.idx

# per-site Fst
realSFS fst print p01.fst.idx > p01.fst

# extracting gene regions out of genome annotations file (gff3)
cat mygenome.gff3 | awk ' $3=="gene"' | cut -f 1,4,5,10 >gene_regions.tab

# extending to plus-minus 2 kb around the gene
awk '{print $1"\t"$2-2000"\t"$3+2000"\t"$4}' gene_regions.tab '> genes.txt

# use fstPerGene.R to compute per-gene Fst

#============= Bayescan: looking for Fst outliers

# Converting vcf (using PGDspider) to Bayescan format: 

# make tab-delimited file called bspops LISTING assignments of individuals (as they are named in the vcf file) to populations, for example:
ind1	pop0
ind2	pop0
ind3	pop1
ind4	pop1

# create a file called vcf2bayescan.spid containing this text:
echo "############
# VCF Parser questions
PARSER_FORMAT=VCF
# Do you want to include a file with population definitions?
VCF_PARSER_POP_QUESTION=true
# Only input following regions (refSeqName:start:end, multiple regions: whitespace separated):
VCF_PARSER_REGION_QUESTION=
# What is the ploidy of the data?
VCF_PARSER_PLOIDY_QUESTION=DIPLOID
# Only output following individuals (ind1, ind2, ind4, ...):
VCF_PARSER_IND_QUESTION=
# Output genotypes as missing if the read depth of a position for the sample is below:
VCF_PARSER_READ_QUESTION=
# Take most likely genotype if "PL" or "GL" is given in the genotype field?
VCF_PARSER_PL_QUESTION=true
# Do you want to exclude loci with only missing data?
VCF_PARSER_EXC_MISSING_LOCI_QUESTION=false
# Select population definition file:
VCF_PARSER_POP_FILE_QUESTION=./bspops
# Only output SNPs with a phred-scaled quality of at least:
VCF_PARSER_QUAL_QUESTION=
# Do you want to include non-polymorphic SNPs?
VCF_PARSER_MONOMORPHIC_QUESTION=false
# Output genotypes as missing if the phred-scale genotype quality is below:
VCF_PARSER_GTQUAL_QUESTION=
# GESTE / BayeScan Writer questions
WRITER_FORMAT=GESTE_BAYE_SCAN
# Specify which data type should be included in the GESTE / BayeScan file  (GESTE / BayeScan can only analyze one data type per file):
GESTE_BAYE_SCAN_WRITER_DATA_TYPE_QUESTION=SNP
############" >vcf2bayescan.spid

# converting vcf (either myresult.vcf from ANGSD or the one from GATK) to bayescan format
java -Xmx1024m -Xms512m -jar ~/bin/PGDSpider_2.0.7.1/PGDSpider2-cli.jar -inputfile myresult.vcf -outputfile Best.bayescan -spid vcf2bayescan.spid 

# launching bayescan (this might take 12-24 hours)
bayescan Best.bayescan -threads=20

# use bayescan_plots.R to examine results

# removing outliers from VCF file
removeBayescanOutliers.pl bayescan=snp.baye_fst.txt vcf=myresult.vcf FDR=0.5 >myresult_nobs.vcf

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

# writing command script with SLURM header (some fields might be different on your cluster, contact your IT people!)
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









