inf=commandArgs(T)

#setwd("~/Dropbox/pipelines2020/auto_ok")
#inf="myresult2.ibsMat"
require(vegan)
ma = as.matrix(read.table(inf))

# reading bams, removing file name trash
bams=scan("bams.nr",what="character") # list of bam files
bams=sub(".bam","",bams)
bams=sub("/.+/","",bams)
dimnames(ma)=list(bams,bams)
# reading sequencing depths (proportion of sites with coverage >5x, for each sample - output of plotQC.R)
qc=read.table("quality.txt",header=F)
names(qc)=c("srr","q")
qc$srr=sub(".bam","",qc$srr)
row.names(qc)=qc$srr
qc=qc[bams,]

# partial PCoA, taking out the variation due to coverage
pp0=capscale(ma~1+Condition(qc$q))

pdf(file=paste(inf,"_pcoa_eigens.pdf",sep=""),height=3.8,width=12)
par(mfrow=c(1,3))
plot(pp0)

# function to calculate how much larder is the dropoff after MDS1+MDS2/2 compared to dropoff in the mid-MDS region
mds1.signal=function(x) {
  N=length(x)
  nulls=round(c(N/2-2,N/2+2),0)
  slope1=(x[1]-x[3])/2
  slope0=(x[nulls[1]]-x[nulls[2]])/4
  return(slope1/slope0)
}

plot(eigsn)

# jackknifing the MDS signal (80% resampling)
mdss=c()
for (i in 1:199) {
  choose=sample(1:nrow(ma),round(nrow(ma)*0.8,0))
  ma1=ma[choose,choose]
  qc1=qc[choose,]
  pp1=capscale(ma1~1+Condition(qc1$q))
  eigsn=pp1$CA$eig/sum(pp1$CA$eig)
  good=mds1.signal(pp1$CA$eig)
  bs=mds1.signal(bstick(length(pp1$CA$eig)))
  # signal of structure = log10 of ratio of real to broken stick MDS1 signals
  mdss=c(mdss,log(good/bs,10))
}
plot(density(mdss),main="log10 ratio of MDS1-2 signal \n to BS model signal, 80% jackknife")
abline(v=median(mdss),col="red")
mtext(paste("median = ",round(median(mdss),2)),side=3,cex=0.8)
dev.off()
save(eigsn,pp1,mdss,file=paste(inf,"_pcoa_eigens.RData",sep=""))
# min, 25%q, median, mean,75%q, max
cat(round(summary(mdss),3))

