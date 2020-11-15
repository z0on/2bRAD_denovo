source('~/Dropbox/pipelines2020/plot_admixture_v5_function.R')
library(RcppCNPy)

# assembling the input table
dir="~/Dropbox/pipelines2020/auto_ok/"
pops="inds2pops" # 2-column tab-delimited table of individual assignments to populations; must be in the same order as samples in the bam list or vcf file.
#------------
setwd(dir)
tbl=data.frame(npyLoad("pcangsd.admix.Q.npy"))
npops=ncol(tbl)
i2p=read.table(paste(dir,pops,sep=""),header=F)

names(i2p)=c("ind","pop")
tbl=cbind(tbl,i2p)
row.names(tbl)=tbl$ind
head(tbl)
tbl$pop=factor(tbl$pop)

ords=plotAdmixture(data=tbl,npops=npops,grouping.method="distance",vshift=0.05,angle=45,cex=0.8)

# recording cluster affiliations
cluster.admix=apply(tbl[,1:npops],1,function(x) {return(paste(which(x>0.25),collapse=".")) })
save(cluster.admix,file=paste(npops,"_clusters.RData",sep=""))
