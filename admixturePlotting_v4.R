# primitive look at admixture data:
# tbl=read.table("rads.2.Q")
# barplot(t(as.matrix(tbl)), col=rainbow(5),xlab="Individual #", ylab="Ancestry", border=NA)

#-----------------
# prettier:

# assembling the input table
dir="" # path to input files
inName="ok_k3.qopt" # name of the input file to plot, output of ngsAdmix or ADMIXTURE run
npops=3
pops="inds2pops" # 2-column tab-delimited table of individual assignments to populations; must be in the same order as samples in the bam list or vcf file.
tbl=read.table(paste(dir,inName,sep=""),header=F)
i2p=read.table(paste(dir,pops,sep=""),header=F)
names(i2p)=c("ind","pop")
tbl=cbind(tbl,i2p)
row.names(tbl)=tbl$ind

head(tbl,20) # this is how the resulting dataset must look

source("plot_admixture_v4_function.R")

# putting populaitons in desired order (edit pop names as needed or skip to plot them alphabetically)
tbl$pop=factor(tbl$pop,levels=c("O","K"))

ords=plotAdmixture(data=tbl,npops=npops,angle=0,vshift=0,hshift=0)
