# ----------------------- INSTALLATIONS --------------------

# --- tools to download SRA datasets from NCBI:

# esearch : see here for installation (on TACC):
https://www.ncbi.nlm.nih.gov/books/NBK179288/ 

# SRA toolkit (on TACC, pick the one for ubuntu 64):
https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit

#--------- cd-hit:

cd
git clone https://github.com/weizhongli/cdhit.git
cd cdhit
make

#---------- Misha's scripts:

cd ~/bin 
# cloning github repositories
git clone https://github.com/z0on/2bRAD_denovo.git
# move scripts to ~/bin from sub-directories
mv 2bRAD_denovo/* . 
# remove now-empty directories
rm -rf 2bRAD_denovo 
cd -

# ------- ANGSD: 

# install xz first from https://tukaani.org/xz/

cd
wget https://tukaani.org/xz/xz-5.2.4.tar.gz --no-check-certificate
tar vxf xz-5.2.4.tar.gz 
cd xz-5.2.4/
./configure --prefix=$HOME/xz-5.2.4/
make
make install

# edit .bashrc:
cd
nano .bashrc
   export LD_LIBRARY_PATH=$HOME/xz-5.2.4/lib:$LD_LIBRARY_PATH
   export LIBRARY_PATH=$HOME/xz-5.2.4/lib:$LIBRARY_PATH
   export C_INCLUDE_PATH=$HOME/xz-5.2.4/include:$C_INCLUDE_PATH
logout
# re-login

# now, install htslib:
cd
git clone https://github.com/samtools/htslib.git
cd htslib
make CFLAGS=" -g -Wall -O2 -D_GNU_SOURCE -I$HOME/xz-5.2.4/include"

# install ANGSD
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

# ---- PCAngsd

# activate your favorite conda environment, then
conda install -c anaconda python
conda install -c anaconda numpy
conda install -c anaconda scipy
conda install -c anaconda cython
cd
git clone https://github.com/Rosemeis/pcangsd.git
cd pcangsd
python setup.py build_ext --inplace
pip3 install -e .

#------------------------------------------------------------------------
#---- CHUNK 1:  getting data

# PRJNA511386: Olympia oysters
# PRJNA430897 : nymphon
# PRJNA343959 : porpoise

export BioProject=PRJNA490084
$HOME/edirect/esearch -db sra -query $BioProject | efetch --format runinfo |cut -d "," -f 1 | grep SRR > $BioProject.SRR && $HOME/edirect/esearch -db sra -query $BioProject | efetch --format runinfo > $BioProject.fullMeta.csv
#esearch -db sra -query $BioProject | efetch --format runinfo | cut -f 1,30 -d "," | grep SRR > $BioProject.srr2sample.csv

>gets
for A in `cat $BioProject.SRR`;do 
echo "fastq-dump-orig.2.10.5 $A">>gets;
done
s2_launcher_creator.py -j gets -n gets -a tagmap -e matz@utexas.edu -t 12:00:00 -w 24 -q normal
getsjob=$(sbatch gets.slurm | grep "Submitted batch job" | perl -pe 's/\D//g')

Q= -log10(Perror)*10
Perror=1% = 0.01 = 1e-2 = 10^(-2)  --> Q=20
P= 0.1% --> 0.001  Q=30


#---- CHUNK 2: processing fastq files (fastq -> bams)

module load cutadapt
module load jellyfish
# module load cd-hit
module load samtools
module load bowtie

export TagLen=100  # TagLen depends on what kind of data you have. 100 is good for most RADs; for 2bRAD, change to 36.
export MatchFrac=0.95 # liberally assuming up to 5% divergence 
export GENOME_REF=all_cc.fasta 
# if you have a real reference genome, put its filename above with full path, don't forget to index it (make a job out of it tho): 
# bowtie2-build mygenome.fasta mygenome.fasta && samtools faidx mygenome.fasta
# you will then skip lines 100-117

# trimming and subsampling to 3M filtered reads max
>trim
for file in *.fastq; do
echo "cutadapt --format fastq -q 15,15 -a AGATCGGA  -m $TagLen -l $TagLen -o ${file/.fastq/}.trim0 $file > ${file}_trimlog.txt && head -12000000 ${file/.fastq/}.trim0 > ${file/.fastq/}.trim && rm ${file/.fastq/}.trim0" >> trim;
done
s2_launcher_creator.py -j trim -n trim -a tagmap -e matz@utexas.edu -t 1:00:00 -w 48 -q normal
trimjob=$(sbatch trim.slurm | grep "Submitted batch job" | perl -pe 's/\D//g')

# converting fastq to fasta (for kmer analysis)
>f2f
for F in `ls *fastq`; do echo "paste - - - - < ${F/.fastq/}.trim | cut -f 1,2 | sed 's/^@/>/' | tr \"\t\" \"\n\" > ${F/.fastq/}.fasta" >>f2f ;done
s2_launcher_creator.py -j f2f -n f2f -a tagmap -e matz@utexas.edu -t 0:05:00 -w 48 -q normal
f2fjob=$(sbatch --dependency=afterok:$trimjob f2f.slurm | grep "Submitted batch job" | perl -pe 's/\D//g')

# analyzing kmers, removing singletons from each file (kmerer.pl)
>jelly
for F in `ls *fastq`; do echo "jellyfish count -m $TagLen -s 100M -t 1 -C ${F/.fastq/}.fasta -o ${F/.fastq/}.jf && jellyfish dump ${F/.fastq/}.jf | kmerer.pl - 2 > ${F/.fastq/}.kmers.fa">>jelly;done
s2_launcher_creator.py -j jelly -n jelly -a tagmap -e matz@utexas.edu -t 0:10:00 -w 24 -q normal
jellyjob=$(sbatch --dependency=afterok:$f2fjob jelly.slurm | grep "Submitted batch job" | perl -pe 's/\D//g')

# merging kmer lists and filtering kmers aiming to get major alleles
# clustering, concatenating into fake genome (10 chromosomes), and indexing it
echo "ls *kmers.fa > all.kf && mergeKmers.pl all.kf minDP=10 minInd=5 > kmers.tab && cat kmers.tab | shuf | awk '{print \">\"\$1\"\n\"\$2}' > all.fasta && cd-hit-est -i all.fasta -o all.clust -aL 1 -aS 1 -g 1 -c $MatchFrac -M 0 -T 0 && concatFasta.pl fasta=all.clust num=10 && bowtie2-build all_cc.fasta all_cc.fasta && samtools faidx all_cc.fasta" >mk
s2_launcher_creator.py -j mk -n mk -a tagmap -e matz@utexas.edu -t 2:00:00 -w 1 -q normal
mergejob=$(sbatch --dependency=afterok:$jellyjob  mk.slurm | grep "Submitted batch job" | perl -pe 's/\D//g')
#mergejob=$(sbatch mk.slurm | grep "Submitted batch job" | perl -pe 's/\D//g')

# mapping, converting to bams, indexing
>maps
for F in `ls *.fastq`; do
REF=all_cc.fasta
echo "bowtie2 --no-unal -x $REF -U ${F/.fastq/}.trim -S ${F/.fastq/}.sam && samtools sort -O bam -o ${F/.fastq/}.bam ${F/.fastq/}.sam && samtools index ${F/.fastq/}.bam">>maps
done
s2_launcher_creator.py -j maps -n maps -a tagmap -e matz@utexas.edu -t 2:00:00 -w 24 -q normal
mapsjob=$(sbatch --dependency=afterok:$mergejob  maps.slurm | grep "Submitted batch job" | perl -pe 's/\D//g')

# quality assessment, removing bams with log(coverage)<3SD
# also imposng minimum number of individuals(MI) a locus must be seen in (genotyping rate cutoff - 50%)
export MinIndPerc=0.5
FILTERSQ='-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minInd $MI'
TODOQ="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2"
echo 'export NIND=`cat bams | wc -l`; export MI=`echo "($NIND*$MinIndPerc+0.5)/1" | bc`' >calc1
echo "ls *.bam > bams && source calc1 && angsd -b bams -r chr1:1-1000000 -GL 1 $FILTERSQ $TODOQ -P 12 -out dd && Rscript ~/bin/plotQC.R prefix=dd">a0
s2_launcher_creator.py -j a0 -n a0 -a tagmap -e matz@utexas.edu -t 2:00:00 -w 1 -q normal
qjob=$(sbatch --dependency=afterok:$mapsjob a0.slurm | grep "Submitted batch job" | perl -pe 's/\D//g')

#------------ examine dd.pdf, decide on MinIndPerc and whether bams.qc is reasonable

# ----------- CHUNK 3: population structure based on common variants

export GENOME_REF=all_cc.fasta
export MinIndPerc=0.8

# initial IBS production, detecting and removing clones (see hctree.pdf and resulting bams.nr)
FILTERS0='-minInd $MI -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -snp_pval 1e-5 -minMaf 0.05 -dosnpstat 1 -doHWE 1 -maxHetFreq 0.5 -hetbias_pval 1e-3 -skipTriallelic 1'
TODO0='-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doPost 1 -doGlf 2'
echo 'export NIND=`cat bams.qc | wc -l`; export MI=`echo "($NIND*$MinIndPerc+0.5)/1" | bc`' >calc1
echo "source calc1 && angsd -b bams.qc -GL 1 $FILTERS0 $TODO0 -P 12 -out myresult && Rscript ~/bin/detect_clones.R bams.qc myresult.ibsMat 0.15">a1
s2_launcher_creator.py -j a1 -n a1 -a tagmap -e matz@utexas.edu -t 2:00:00 -w 1 -q normal
# a1job=$(sbatch --dependency=afterok:$qjob a1.slurm | grep "Submitted batch job" | perl -pe 's/\D//g')
a1job=$(sbatch a1.slurm | grep "Submitted batch job" | perl -pe 's/\D//g')

# if "highly similar samples" were reported in a1.e* file, examine hctree.pdf and possibly rerun
# Rscript ~/bin/detect_clones.R bams.qc myresult.ibsMat 0.15
# with higher or lower cutoff instead of 0.15

# final IBS production, assessing "PCA structure" (excess of MDS1-2 signal relative to broken stick model)
# minMaf filter set to 3/2N (minimal AF for an allele to be guaranteed to be found in >1 individual, to minimize ascertainment bias)
FILTERS1='-minInd $MI2 -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -snp_pval 1e-5 -minMaf $MAF -dosnpstat 1 -doHWE 1 -maxHetFreq 0.5 -hetbias_pval 1e-3 -skipTriallelic 1'
TODO1='-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doPost 1 -doGlf 2'
echo 'cat bams.nr | sort > bams.NR && mv bams.NR bams.nr && export NIND2=`cat bams.nr | wc -l`; export MI2=`echo "($NIND2*$MinIndPerc+0.5)/1" | bc`; export MAF=`echo "3/(2*$NIND2)" | bc -l`' >calc2
echo "source calc2 && angsd -b bams.nr -GL SeqTools1 $FILTERS1 $TODO1 -P 12 -out myresult2 && Rscript ~/bin/pcaStructure.R myresult2.ibsMat">a2
s2_launcher_creator.py -j a2 -n a2 -a tagmap -e matz@utexas.edu -t 2:00:00 -w 1 
-q normal
a2job=$(sbatch --dependency=afterok:$a1job a2.slurm | grep "Submitted batch job" | perl -pe 's/\D//g')
#a2job=$(sbatch  a2.slurm | grep "Submitted batch job" | perl -pe 's/\D//g')

# PCAngsd: admixture, SNP covariance, kinship, inbreeding (based on high-frequency SNPs! better measure is individual heterozygosity, l 206-209)

module load python2
echo 'python ~/pcangsd/pcangsd.py -beagle myresult2.beagle.gz -admix -o pcangsd -inbreed 2 -kinship -selection -threads 12' >adm
s2_launcher_creator.py -j adm -n adm -a tagmap -e matz@utexas.edu -t 0:10:00 -w 1 -q normal
admjob=$(sbatch --dependency=afterok:$a2job adm.slurm | grep "Submitted batch job" | perl -pe 's/\D//g')

# ------- genetic diversity stats based on all well-genotyped sites (including invariants)
# (note the argument -r chr5 on line 199, ie only using chr 5 - we assume that this region is enough to reresent the whole genome)

export GENOME_REF=all_cc.fasta
export MinIndPerc=0.8

FILTERS='-minInd $MI2 -uniqueOnly 1 -skipTriallelic 1 -minMapQ 20 -minQ 20 -doHWE 1 -maxHetFreq 0.5 -hetbias_pval 1e-3'
TODO="-doSaf 1 -anc $REF -ref $GENOME_REF -doMajorMinor 1 -doMaf 1 -dosnpstat 1 -doPost 1 -doGlf 2"
echo 'export NIND2=`cat bams.nr | wc -l`; export MI2=`echo "($NIND2*$MinIndPerc+0.5)/1" | bc`' >calc2
echo "source calc2 && angsd -b bams.nr -r chr5 -GL 1 -P 12 $FILTERS $TODO -out chr5 && realSFS chr5.saf.idx -P 12 -fold 1 > chr5.sfs && realSFS saf2theta chr5.saf.idx -outname chr5 -sfs chr5.sfs -fold 1 && thetaStat do_stat chr5.thetas.idx -outnames chr5 && grep \"chr\" chr5.pestPG | awk '{ print \$4/\$14}' >piPerChrom">sfsj
s2_launcher_creator.py -j sfsj -n sfsj -t 2:00:00 -e matz@utexas.edu -w 1 -a tagmap -q normal
sfsjob=$(sbatch --dependency=afterok:$a1job sfsj.slurm | grep "Submitted batch job" | perl -pe 's/\D//g')
#sfsjob=$(sbatch sfsj.slurm | grep "Submitted batch job" | perl -pe 's/\D//g')


# individual heterozygosities (proportion of heterozygotes across SNPs that pass sfsjob filters)
echo "Rscript ~/bin/heterozygosity_beagle.R chr5.beagle.gz >indHets" >bg
s2_launcher_creator.py -j bg -n bg -a tagmap -e matz@utexas.edu -t 12:00:00 -w 1 -q normal
sbatch --dependency=afterok:$sfsjob bg.slurm
# heterozygosity_beagle.R script (by Nathaniel Pope) outputs *_zygosity.RData R data bundle containing AFS (rows) for each individual (columns). The proportion of heterozygotes is the second row. Individual heterozygosities are also printed out to STDOUT (in the case above, saved to text file indHets)

# relatedness with NgsRelate
echo 'export NIND2=`cat bams.nr | wc -l`; export NS=``zcat g3.mafs.gz | wc -l`' >calc3
echo 'source calc3 && zcat g3.mafs.gz | cut -f5 |sed 1d >freq && ngsRelate  -g g3.glf.gz -n $NIND -f freq >g3.relatedness' >rel
s2_launcher_creator.py -j rel -n rel -a tagmap -e matz@utexas.edu -t 2:00:00 -w 1 -q normal
reljob=$(sbatch rel.slurm | grep "Submitted batch job" | perl -pe 's/\D//g')








#### the following should probably be retired - inbreeding is better done in the form of individual heterozygosities, and relatedness ~ kinship from PCAngsd. (did not explore this yet, tho)


#---------------- if you want to analyze LD, relatedness, or inbreeding :

# rerunning with -doGeno 8 and -doGlf 3 for ngsLD, ngsRelate and ngsF
FILTERS1='-minInd $MI2 -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -snp_pval 1e-5 -minMaf 0.05 -dosnpstat 1 -doHWE 1 -maxHetFreq 0.5 -hetbias_pval 1e-3 -skipTriallelic 1'
TODO1='-doMajorMinor 1 -doMaf 1 -doCounts 1 -doPost 1 -doGlf 3'
echo 'cat bams.nr | sort > bams.NR && mv bams.NR bams.nr && export NIND2=`cat bams.nr | wc -l`; export MI2=`echo "($NIND2*$MinIndPerc+0.5)/1" | bc`' >calc2
echo "source calc2 && angsd -b bams.nr -GL 1 $FILTERS1 $TODO1 -P 12 -out g3">g3
s2_launcher_creator.py -j g3 -n g3 -a tagmap -e matz@utexas.edu -t 2:00:00 -w 1 -q normal
g3job=$(sbatch g3.slurm | grep "Submitted batch job" | perl -pe 's/\D//g')


# individual inbreeding coeffs
echo 'zcat g3.glf.gz | ngsF --glf - --n_ind 44 --n_sites 909052 --out inbr' >nf
s2_launcher_creator.py -j nf -n nf -a tagmap -e matz@utexas.edu -t 2:00:00 -w 1 -q normal
ngsFjob=$(sbatch --dependency=afterok:$g3job nf.slurm | grep "Submitted batch job" | perl -pe 's/\D//g')
#ngsFjob=$(sbatch nf.slurm | grep "Submitted batch job" | perl -pe 's/\D//g')

zcat g3.mafs.gz | cut -f5 |sed 1d >freq

# relatedness with NgsRelate
echo 'export NIND2=`cat bams.nr | wc -l`; export NS=``zcat g3.mafs.gz | wc -l`' >calc3
echo 'source calc3 && zcat g3.mafs.gz | cut -f5 |sed 1d >freq && ngsRelate  -g g3.glf.gz -n $NIND2 -f freq >g3.relatedness' >rel
s2_launcher_creator.py -j rel -n rel -a tagmap -e matz@utexas.edu -t 2:00:00 -w 1 -q normal
reljob=$(sbatch rel.slurm | grep "Submitted batch job" | perl -pe 's/\D//g')



