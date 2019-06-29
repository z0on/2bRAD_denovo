
#setwd("~/Dropbox/mega2019/idotea_2019/ADMIX/")
source("plot_admixture_v4_function.R")

# lat-lon info, for pops sorting
meta=read.csv("~/Dropbox/mega2019/idotea_2019/Idotea environment 171005.csv")
i2p=read.table("~/Dropbox/mega2019/idotea_2019/ADMIX/inds2pops",sep="\t")
names(i2p)=c("ind","pop")
mordi=merge(i2p, meta[,c(1:3,5:8)],by="pop",all.x=T)
head(mordi)
mordi$ind==i2p$ind
byLon=i2p$ind[order(mordi$lon)]

# assembling the input table
dir="~/Dropbox/mega2019/idotea_2019/ADMIX/" # path to input files
inName="idotea600_k6.qopt" # name of the input file to plot, output of ngsAdmix or ADMIXTURE run
pops="inds2pops" # 2-column tab-delimited table of individual assignments to populations; must be in the same order as samples in the bam list or vcf file.

#------------

npops=as.numeric(sub(".+(\\d+)\\..+","\\1",inName))
tbl=read.table(paste(dir,inName,sep=""),header=F)
i2p=read.table(paste(dir,pops,sep=""),header=F)

names(i2p)=c("ind","pop")
tbl=cbind(tbl,i2p)
row.names(tbl)=tbl$ind
head(tbl[byLon,],20) # this is how the resulting dataset must look

# putting populaitons in desired order (edit pop names as needed or skip to plot them alphabetically)
# tbl$pop=factor(tbl$pop,levels=c("O","K"))

ords=plotAdmixture(data=tbl,npops=npops,ord=byLon,grouping.method="distance",hshift=9,vshift=0.05,cex=0.6,srt=45)

# recording cluster affiliations
cluster.admix=apply(tbl[,1:npops],1,function(x) {return(paste(which(x>0.25),collapse=".")) })
#save(cluster.admix,file=paste(inName,"_clusters.RData",sep=""))
