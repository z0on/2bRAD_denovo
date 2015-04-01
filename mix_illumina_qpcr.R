# This script calculates and plots the number of microliters of each sample 
# that should be mixed to achieve equal concentration,
# based on qPCR data using Illumina standard library (truSeq?) primers
# Each sample should be qPCR'ed as two dilutions, a more concentrated one (1/10 or 1/16 of the gel extracted DNA), 
# and another 10-fold down from this one; two replicates of each dilution (replicate the dilutions themselves!).

# The model will evaluate amplification efficiency of each sample and infer its intercept, which is what is needed to 
# calculate relative concentrations

# The intercept is computed as a random effect rather than fixed effect, to avoid occasional bad overshooting due to 
# small number of data points for each sample.

# the result will be written to a tab-delimited file "illumina_mixing_chart.txt", open it in excel

# note that by defaiult the script outputs 75% confidence intervals (I feel this is more practical). If you want 95%, specify prob=95 in a call
# to the function dna.mixing (the second line from the end)

# The script requires requires MCMCglmm package. If you don't have it installed, un-remark the next line and execute it (need to do it just once):
# install.packages("MCMCglmm")

# ready? start executing the chunks between #--------- lines one by one, paying attention to the notes (sometimes you will need to remark some 
# lines and un-remark others, depending on your data structure). 
# To execute a chunk of R code, mark it and hit command-return on Mac, control-enter on PC

# Mikhail V. Matz, UT Austin, matz@utexas.edu
# August 30, 2012

# ------------------
# press command-D on Mac or Alt-F-C on a PC and navigate to the directory containing your data file and scripts
#-------------------
dat=read.csv(file.choose())  # open your data file, or illumina_mix_data.csv for practicing and formatting example
head(dat)
# there must be four columns: sam (sample name), lane (intended HiSeq lane), conc (DNA dilution; use 0.0625 for 1/16 and 0.00625 for 1/160), and ct (qPCR result for this sam-conc combination). There must be at least two technical replicates for each combination of sam-conc (i.e. two rows with teh same sam and conc and different ct values). If all samples are to be mixed in the same lane, enter '1' throughout the lane column.
# the order of columns or rows does not matter, the names of the columns do matter (and they are case sensitive)

dat$lane=as.factor(dat$lane)
dat$logc=log(dat$conc,2)

#-------------------

source("dna.mixing.R")
library(MCMCglmm)

# if not all samples have at least two replicates for each dilution, un-remark the next line and remark the next-next one
#lmf1=MCMCglmm(ct~logc,random=~sam,data=dat,verbose=F,pr=T)
lmf1=MCMCglmm(ct~logc,random=~idh(1+logc):sam,data=dat,pr=T)
slop=mean(lmf1$Sol[,2])
blups=colMeans(lmf1$Sol)
slops=blups[grep("sam.logc.sam",names(blups))]
names(slops)=sub("sam.logc.sam.","",names(slops))
inters1=blups[grep("sam.\\(Intercept\\).sam",names(blups))]

es=2^(-1/(slop+slops))
hist(es)
# this the histogram of the amplification efficiencies. Ideally they all should be within 1.9-2.0 range; but sometimes things get weird (efficiencies 1.7 or 2.5) for unknown reason. Fortunately, this does not seem to affect the mixing, as long as all the samples show similar efficiencies (i.e., within 0.1 of each other) 

#----------------------

par(mfrow=c(floor(sqrt(length(levels(dat$lane)))),ceiling(sqrt(length(levels(dat$lane))))+1))
result=dna.mixing(lmf1,dat) # this plot can be scaled by dragging the lower-right corner, and saved as pdf when pretty enough (click on the plot window, then click File-SaveAs in the menu)
result$lane=as.numeric(result$lane)
result=round(result[,1:4],1)
result
# look at the first column only, ul2mix - this gives the number of microliters of the original eluate to put into the lane-specific mix to get equal representation of all samples within the lane. Samples of very low concentration might require more than there is of the eluate (100 ul or more); in that case, consider diluting the most-concentrated samples from the same lane (i.e., take the equivalent of 0.1 ul instead of 1 ul of the most concentrated one, which will lower the requirement 10-fold across the whole lane), remaking the low-concentrated ones from scratch, or just dropping them. 
write.table(result, file="illumina_mixing_chart.txt",quote=F,sep="\t") # this table can be opened in Excel

