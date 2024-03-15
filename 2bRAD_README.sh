2bRAD de novo AND reference-based walkthrough
July 2021
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

------- cutadapt: 

cdh
pip install --user cutadapt
cp .local/bin/cutadapt ~/bin

------- cd-hit:

git clone https://github.com/weizhongli/cdhit.git
cd cd-hit
make

------- Moments: 

cd
git clone https://bitbucket.org/simongravel/moments.git 
cd moments
python setup.py build_ext --inplace

# add this to .bashrc, section 2:
  export PYTHONPATH=$PYTHONPATH:$HOME/moments
# re-login

# NOTA BENE: on ls5, must replace 'set_axis_bgcolor' with 'set_facecolor' in moments script ModelPlot.py (in moments directory, typically, $HOME/moments/moments/)

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

-------  ngsRelate :
cd 
git clone https://github.com/ANGSD/NgsRelate.git
cd NgsRelate
make HTSSRC=../htslib
cp ngs* ~/bin/
cd

------- ngsF  : (inbreeding coefficients)

git clone https://github.com/fgvieira/ngsF.git
module load gsl
cd ngsF

nano Makefile
add -I${TACC_GSL_INC}  to CC and CXX macros (CFLAGS= ...);
and -L${TACC_GSL_LIB} to the 'LIB = ...' line.

export PKG_CONFIG_PATH=/opt/apps/intel18/gsl/2.2.1/lib/pkgconfig/
make HTSSRC=../htslib

------ ngsLD :

cd 
git clone https://github.com/fgvieira/ngsLD.git
cd ngsLD

nano Makefile
add -I${TACC_GSL_INC}  to CC and CXX macros (CFLAGS= ...);
and -L${TACC_GSL_LIB} to the 'LIB = ...' line.

module load gsl
export PKG_CONFIG_PATH=/opt/apps/intel18/gsl/2.2.1/lib/pkgconfig/
make
cp ngsLD ~/bin

-------  stairwayPlot :

# project page: https://sites.google.com/site/jpopgen/stairway-plot
cdw
# get version from June 2016 (v2beta2)
wget https://www.dropbox.com/s/toxnlvk8rhe1p5h/stairway_plot_v2beta2.zip
unzip stairway_plot_v2beta2.zip
mv stairway_plot_v2beta2 stairway_plot_v2beta

# ---- PCAngsd

cd
module load python2
git clone https://github.com/Rosemeis/pcangsd.git
cd pcangsd/
python setup.py build_ext --inplace

-------  ADMIXTURE

cd ~/bin/
wget http://software.genetics.ucla.edu/admixture/binaries/admixture_linux-1.3.0.tar.gz --no-check-certificate
tar vxf admixture_linux-1.3.0.tar.gz 
mv admixture_linux-1.3.0/admixture .
cd -

-------  plink 1.9:
cd ~/bin
wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20190617.zip
unzip plink_linux_x86_64_20190617.zip
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
#git clone https://github.com/z0on/2bRAD_GATK.git
# move scripts to ~/bin from sub-directories
#mv 2bRAD_GATK/* . 
mv 2bRAD_denovo/* . 
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

# skip this unless you plan to use GATK (hard-call genotyping, with >10x coverage)
export GENOME_DICT=$WORK/db/mygenome.dict 
java -jar $TACC_PICARD_T_DIR/picard.jar CreateSequenceDictionary R=$GENOME_FASTA  O=$GENOME_DICT

#==================
# Step 1: Splitting by in-read barcode, deduplicating and quality-filtering the reads

# creating a file of commands to run (assuming reads are in fastq files, one file per sample.
# (if samples were spread across several lanes, concatenate them first using ngs_concat.pl) 
2bRAD_trim_launch_dedup.pl fastq > trims

# Note: use this command instead of the one above if you have 2bRAD libraries made with oligos ordered before March 2021:
2bRAD_trim_launch_dedup_old.pl fastq > trims
# use this command instead of the one above if you have 2bRAD libraries without degenerate tag but with 4-base in-line barcode:
2bRAD_trim_launch.pl fastq barcode2=4 > trims
# And if you have the original-style libraries without degenerate tag and without in-line barcode, use this:
2bRAD_trim_launch.pl fastq > trims
# AND IF you sequenced your double-barcoded, deduplicatable libraries on HiSeq 4000 alone
# (resulting in poor quality at restriction site and adaptor bases) and you have used BcgI enzyme, use this:
2bRAD_trim_launch_dedup2.pl fastq > trims
# Once on Nextseq we saw annoying Ns in the upstream adaptor (in place of invariant CC there), if you encounter that, use this:
2bRAD_trim_launch_dedup_N2.pl fastq > trims


# the file trims now contains a long list of commands, one per fastq file, that we will need to execute. TACC uses Launcher module to enable parallele launch of multiple commands: https://github.com/TACC/launcher 
# I very highly recommend it - it should work on any SLURM-based system. Otherwise, if you don't know how to run hundreds of parallel jobs on your cluster, consult your IT support - there must be a way. 
# If push comes to shove, you can always execute all these commands serially (one by one), by saying 
bash trims

# do we have expected number of *.tr0 files created?
ls -l *.tr0 | wc -l

# quality filtering using cutadapt (see installation above)
module load cutadapt 

# for de novo analysis: removing reads with qualities at ends less than Q15
>trimse
for file in *.tr0; do
echo "cutadapt --format fastq -q 15,15 -m 36 -o ${file/.tr0/}.trim $file > ${file}_trimlog.txt" >> trimse;
done

# for reference-based analysis: trimming poor quality bases off ends:
>trimse
for file in *.tr0; do
echo "cutadapt --format fastq -q 15,15 -m 25 -o ${file/.tr0/}.trim $file > ${file}_trimlog.txt" >> trimse;
done

# execute all commands in trimse file (serial or parallel using Launcher, if your system allows) 

# do we have expected number of *.trim files created?
ls -l *.trim | wc -l

#====================
# denovo RAD business (skip to next "#========" if you have genome reference):

# 'uniquing' ('stacking') individual trimmed fastq reads:
ls *.trim | perl -pe 's/^(.+)$/uniquerOne.pl $1 >$1\.uni/' >unii

# execute all commands written to unii...

# Done! do you have .uni for all your samples?... 
ls -l *.uni | wc -l  

# collecting common tags (= major alleles)
# merging uniqued files (set minInd to >10, or >10% of total number of samples, whichever is greater)
mergeUniq.pl uni minInd=10 >all.uniq

# discarding tags that have more than 7 observations without reverse-complement
awk '!($3>7 && $4==0) && $2!="seq"' all.uniq >all.tab

# creating fasta file out of merged and filtered tags:
awk '{print ">"$1"\n"$2}' all.tab > all.fasta

# clustering allowing for up to 3 mismatches (-c 0.91); the most abundant sequence becomes reference
cd-hit-est -i all.fasta -o cdh_alltags.fas -aL 1 -aS 1 -g 1 -c 0.91 -M 0 -T 0  

#------------
# making fake reference genome (of 30 chromosomes) out of major-allele tags
# need bowtie2 and samtools for indexing

concatFasta.pl fasta=cdh_alltags.fas num=20

# formatting fake genome
export GENOME_FASTA=cdh_alltags_cc.fasta
export GENOME_DICT=cdh_alltags_cc.dict 
bowtie2-build $GENOME_FASTA $GENOME_FASTA
samtools faidx $GENOME_FASTA

# skip this one if you don't plan to use GATK:
java -jar $WHERE_PICARD_IS/picard.jar CreateSequenceDictionary R=$GENOME_FASTA  O=$GENOME_DICT

#==============
# Mapping reads to reference (reads-derived fake one, or real) and formatting bam files 

# for denovo: map reads to fake genome: 
GENOME_FASTA=cdh_alltags_cc.fasta

# for reference-based: 
GENOME_FASTA=mygenome.fasta

# mapping with --local option, enables clipping of mismatching ends (guards against deletions near ends of RAD tags)
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
module load samtools
>s2b
for file in *.sam; do
echo "samtools sort -O bam -o ${file/.sam/}.bam $file && samtools index ${file/.sam/}.bam">>s2b;
done

# run all commands listed in s2b file
ls *bam | wc -l  # should be the same number as number of trim files

# BAM files are the input into various genotype calling / popgen programs, this is the main interim result of the analysis. Archive them.

#===================== A  N  G  S  D =====================

# "FUZZY genotyping" with ANGSD - without calling actual genotypes but working with genotype likelihoods at each SNP. Optimal for low-coverage data (<10x).
# if your coverage is >10x, go to GATK section below

# install ANGSD first (see Installations section above)

# listing all bam filenames 
ls *bam >bams

#----------- assessing base qualities and coverage depth

# angsd settings:
# -minMapQ 20 : only highly unique mappings (prob of erroneous mapping =< 1%)
# -baq 1 : realign around indels (not terribly relevant for 2bRAD reads mapped with --local option) 
# -maxDepth : highest total depth (sum over all samples) to assess; set to 10x number of samples
# -minInd : the minimal number of individuals the site must be genotyped in. Reset to 50% of total N at this stage.

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -maxDepth 1000 -minInd 1000"

# T O   D O : 
TODO="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2"

# in the following line, -r argument is one chromosome or contig to work with (no need to do this for whole genome as long as the chosen chromosome or contig is long enough, ~1 Mb. When mapping to a real genome, consider chr1:1-1000000 )
# (look up lengths of your contigs in the header of *.sam files)
angsd -b bams -r chr1 -GL 1 $FILTERS $TODO -P 1 -out dd 

# summarizing results (using modified script by Matteo Fumagalli)
Rscript ~/bin/plotQC.R prefix=dd
# proportion of sites covered at >5x:
cat quality.txt

# scp dd.pdf to laptop to look at distribution of base quality scores, fraction of sites in each sample passing coverage thresholds, and fraction of sites passing genotyping rates cutoffs. Use these to guide choices of -minQ,  -minIndDepth and -minInd filters in subsequent ANGSD runs

#--------------- population structure

# Note: PCA and Admixture are not supposed to be run on data that contain clones or genotyping replicates. For PCA, these can be removed without rerunning ANGSD from the IBS distance matrix; but for ngsAdmix ANGSD must be rerun.

# Generating genotype likelihoods from highly confident (non-sequencing-error) SNPs
# set minInd to 75-80% of your total number of bams
# if you expect very highly differentiated populations with nearly fixed alternative alleles, remove '-hwe_pval 1e-5' form FILTERS
# -doGeno 8 : genotype likelihood format setting for ngsLD; if you want to run PCA, use -doGeno 32 (but I recommend using ibsMat for all ordination work)
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 1000 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBCF 1 -doPost 1 -doGlf 2"

# Starting angsd with -P the number of parallel processes. Funny but in many cases angsd runs faster on -P 1
angsd -b bams -GL 1 $FILTERS $TODO -P 1 -out myresult

# how many SNPs?
NSITES=`zcat myresult.mafs.gz | wc -l`
echo $NSITES

# ------- PCAngsd: estimating admixture, kinship, and SNP covariance (note: inbreeding is better estimated as individual heterozygosities later, based on all sites, to avoid ascertainment bias)

module load python2
python ~/pcangsd/pcangsd.py -beagle myresult.beagle.gz -admix -o pcangsd -inbreed 2 -kinship -selection -threads 12

# scp bams, myresult.ibsMat, and pcangsd* files to laptop, plot it in R:
# use admixturePlotting_pcangsd.R to plot (will require a table of bams vs sampling locations)

# scp *Mat, *qopt and bams files to laptop, use angsd_ibs_pca.R to plot PCA and admixturePlotting_v4.R to plot ADMIXTURE

# ----- individual heterozygosity by running angsd -dosaf on individual bams
# using all sites (variable and invariable)

# for denovo: fake genome: 
GENOME_FASTA=cdh_alltags_cc.fasta
# for reference-based: 
GENOME_FASTA=mygenome.fasta

# adjust setMinDepthInd as desired. Currently set to only look at sites with x10+ coverage.
FILTERS='-minInd 1 -setMinDepthInd 10 -uniqueOnly 1 -minMapQ 30 -minQ 30'
>hets
>mybams.het
for F in `cat bams`; do
echo "angsd -i $F -anc $GENOME_FASTA $FILTERS -GL 1 -dosaf 1 -out ${F/.bam/} && realSFS ${F/.bam/}.saf.idx >${F/.bam/}.ml | awk -v file=$F '{print file\"\t\"(\$1+\$2+\$3)\"\t\"\$2/(\$1+\$2+\$3)}' ${F/.bam/}.ml >>mybams.het">>hets;
done
# execute all commands in hets; it will produce a tab-delimited table mybams.het:
# [bam filename]   [total number of sites passing filters]   [heterozygosity]

# ------ LD: (use rEM column for LD networks method, to look for signatures of polygenic selection):
NS=`zcat myresult.geno.gz | wc -l`
NB=`cat bams | wc -l`
zcat myresult.mafs.gz | tail -n +2 | cut -f 1,2 > mc1.sites
module load gsl
ngsLD --geno myresult.geno.gz --probs 1 --n_ind $NB --n_sites $NS --max_kb_dist 0 --pos mc1.sites --out myresult.LD --n_threads 12 --extend_out 1


# relatedness (must run ANGSD with option '-doGlf 3' to make this work)
# (column "rab" in the result is relatedness coefficient, Fa and Fb are individual inbreeding coefficients - beware of ascentainment bias tho!):
zcat myresult.mafs.gz | cut -f5 |sed 1d >freq
NIND=`cat bams | wc -l`
ngsRelate -f freq -g myresult.glf.gz -n $NIND -z bams >relatedness


#======= ANDSD on all well-genotyped sites: genetic diversity stats, SFS for demographic analysis

# make separate files listing bams for each population (without clones and replicates)
# assume we have two populations, pop0 and pop1, 20 individuals each, with corresponding bams listed in pop0.bams and pop1.bams

# NOTE: for speed and memory sake, it might make sense to run the following on a subset of the data,
# specifying -r or -rf argument for angsd to look only at a fraction of the reference genome.

# generating list of filtered SNP sites for SFS production (note: no filters that distort allele frequency!):
# sb - strand bias filter; only use for 2bRAD, GBS or WGS (not for ddRAD or RADseq)
# hetbias - detects weird heterozygotes because they have unequal representation of alleles
# maxHetFreq - filters out lumped paralogs 
# set minInd to 80% of all your individuals (depending on the results from quality control step)
# If the goal is genome-wide diversity statistics, consider running this on a fraction of genome (a few Mb) - use angsd options -r or -rf

# ---- stage 1: finding well-genotytped sites (both variable and invariable) that are shared between two populations.
# we will run filtering command on each population, WITHOUT allele-frequency distorting filters (-minMaf and -snp_pval)

# assuming we have text files `pop0.bams` and `ppop1.bams` listing bam files for two populations:
# computing setting for the -minInd filter:
export GenRate=0.75 # desired genotyping rate
export N1=`wc -l pop0.bams | cut -f 1 -d " "`
export N2=`wc -l pop1.bams | cut -f 1 -d " "`
export MI1=`echo "($N1*$GenRate+0.5)/1" | bc`
export MI2=`echo "($N2*$GenRate+0.5)/1" | bc`

FILTERS='-uniqueOnly 1 -skipTriallelic 1 -minMapQ 30 -minQ 30 -maxHetFreq 0.5 -hetbias_pval 1e-3'
# add `-sb_pval 1e-3` (strand bias) to FILTERS if you have 2bRAD, GBS, or WGS data. Other types of RAD only sequence one strand so -sb_pval filter would remove everything.
GENOME_REF=mygenome.fasta # reference to which the reads were mapped
TODO="-doHWE 1 -doSaf 1 -doMajorMinor 1 -doMaf 1 -doPost 2 -dosnpstat 1 -doGeno 11 -doGlf 2 -anc $GENOME_REF -ref $GENOME_REF"
angsd -b pop0.bams -GL 1 -p 4 -minInd $MI1 $FILTERS $TODO -underFlowProtect 1 -out pop0
angsd -b pop1.bams -GL 1 -P 4 -minInd $MI2 $FILTERS $TODO -underFlowProtect 1 -out pop1

# collecting and indexing filter-passing sites
zcat pop0.mafs.gz | cut -f 1,2 | tail -n +2 | sort >pop0.sites
zcat pop1.mafs.gz | cut -f 1,2 | tail -n +2 | sort >pop1.sites

# collecting and indexing common sites:
comm -12 pop0.sites pop1.sites | sort -V >allSites
angsd sites index allSites

# listing "regions"
cat allSites | cut -f 1 | uniq >regions

# rerun saf generation for each pop using common sites
GENOME_REF=mygenome.fasta
angsd -b pop0.bams -rf regions -sites allSites -GL 1 -P 8 -doSaf 1 -anc $GENOME_REF -underFlowProtect 1 -out pop0s
angsd -b pop1.bams -rf regions -sites allSites -GL 1 -P 8 -doSaf 1 -anc $GENOME_REF -underFlowProtect 1 -out pop1s


#===================== 2d AFS analysis using Moments

# This part has it own github page now: https://github.com/z0on/AFS-analysis-with-moments
# If you just want to generate nice "bagged" AFS for plotting:

realSFS pop0s.saf.idx pop1s.saf.idx -ref $GENOME_REF -anc $GENOME_REF -bootstrap 5 -P 1 -resample_chr 1 >p12

# computing SFS dimensions
export N1=`wc -l pop0.bams | cut -f 1 -d " "`
export N2=`wc -l pop1.bams | cut -f 1 -d " "`
export NG1=`echo "($N1*2)+1" | bc`
export NG2=`echo "($N2*2)+1" | bc`

# averaging bootstraps, writing sfs file with header
echo "$NG1 $NG2">p12.sfs
cat p12 | awk '{for (i=1;i<=NF;i++){a[i]+=$i;}} END {for (i=1;i<=NF;i++){printf "%.3f", a[i]/NR; printf "\t"};printf "\n"}' >> p12.sfs

# To plot and project this SFS to fewer samples, use python script 2dAFS.py from this location: 
# https://github.com/z0on/AFS-analysis-with-moments/tree/master/multimodel_inference/py3_v1
# Requires matplotlib, pylab, numpy, sys, and moments packages

# This SFS can also be plotted in R simply, like this: 
'''R
sfs2matrix=function(sfs,n1,n2,zero.ends=TRUE) {
  dd=matrix(ncol=2*n1+1,nrow=2*n2+1,sfs)
  if(zero.ends==TRUE) { dd[1,1]=dd[2*n2+1,2*n1+1]=0 }
  return(apply(dd,2,rev))
}
s=scan("p12.sfs")
n1=(s[1]-1)/2
n2=(s[2]-1)/2
ssfs=sfs2matrix(s[-c(1:2)],n1,n2)
plot(raster(log(ssfs+0.1,10)))
'''

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

# Converting vcf or bcf (using PGDspider) to Bayescan format: 

# make tab-delimited file called bspops LISTING assignments of individuals (as they are named in the vcf file) to populations, for example:
ind1	pop0
ind2	pop0
ind3	pop1
ind4	pop1


# create a file called vcf2bayescan.spid containing this text:
echo "############
# VCF Parser questions (replace VCF with BCF in the line below if you used ANGSD to make bcf instead of vcf)
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
VCF_PARSER_EXC_MISSING_LOCI_QUESTION=true
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
java -Xmx1024m -Xms512m -jar ~/bin/PGDSpider_2.0.7.1/PGDSpider2-cli.jar -inputfile OKbs.vcf -outputfile Best.bayescan -spid vcf2bayescan.spid 

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









