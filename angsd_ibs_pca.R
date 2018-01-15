setwd("~/Dropbox/Documents/ecogeno2018/")
bams=read.table("bams")[,1] # list of bam files
goods=c(1:length(bams))

# reading table of pairs of replicates (tab-delimited) - skip if there are no clones
clonepairs=read.table("clonepairs.tab",sep="\t")
repsa= clonepairs[,1]
repsb= clonepairs[,2]
# removing "b" replicates
goods=which(!(bams %in% repsb))

#--------------------
# loading individual to population correspondences
i2p=read.table("inds2pops",sep="\t") # 2-column tab-delimited table of individual assignments to populations; must be in the same order as samples in the bam list or vcf file.
row.names(i2p)=i2p[,1]
i2p=i2p[goods,]
site=i2p[,2]

# settign up colors for plotting
palette(rainbow(length(unique(site))))
colors=as.numeric(as.factor(site))
colpops=as.numeric(as.factor(sort(unique(site))))
#--------------------
# covariance / PCA 

library(vegan)
# choose either of the following two covarince matrices:
co = as.matrix(read.table("ok.covMat")) # covariance based on single-read resampling
#co = as.matrix(read.table("ok.covar")) # covariance by ngsCovar
co =co[goods,goods]
dimnames(co)=list(bams[goods],bams[goods])

# PCoA and CAP (constranied analysis of proximities)  
conds=data.frame(cbind(site))
pp0=capscale(as.dist(1-cov2cor(co))~1) # PCoA
pp=capscale(as.dist(1-cov2cor(co))~site,conds) # CAP

# significance of by-site divergence, based on 1-correlation as distance
adonis(as.dist(1-cov2cor(co))~site,conds)

# eigenvectors: how many are interesting?
plot(pp0$CA$eig) 

axes2plot=c(1,2)  
quartz()
cc=pp0
plot(cc,choices=axes2plot,type="n") # choices - axes to display
points(cc,choices=axes2plot,pch=19,col=colors)
#ordihull(cmd,choices= axes2plot,groups= conds$grp,draw="polygon",col=1+as.numeric(unique(as.factor(conds$grp))),label=T)
ordispider(cc,choices= axes2plot,groups=conds$site,col="grey80")
ordiellipse(cc,choices= axes2plot,groups= conds$site,draw="polygon",col=colpops,label=T)

# unscaled, to identify outliers
n2identify=2
cmd=pp0
# unscaled, to identify outliers
plot(cmd$CA$u[,axes2plot],pch=19,col=colors)
ordispider(cmd$CA$u[,axes2plot],groups=conds$site,col="grey80")
ordiellipse(cmd$CA$u[,axes2plot],groups= conds$site,draw="polygon",col=colpops,label=T)
identify(cmd$CA$u[,axes2plot],labels=colnames(co),n=n2identify,cex=0.7)
#-------------
# t-SNE:  machine learning to identify groups of samples 
# based on genotypes' correlations

library(Rtsne)
library(vegan)
library(adegenet)
quartz()

# perplexity:  expected number fo neighbors. Set to 0.5x N(samples per pop)
perp=15
rt = Rtsne(as.dist(1-cov2cor(co)), perplexity=perp,max_iter=2,is_distance=T)
for (i in 1:250){
	rt = Rtsne(as.dist(1-cov2cor(co)), perplexity=perp,max_iter=10,Y_init=rt$Y,is_distance=T)
	plot(rt$Y,col=colors,pch=16,cex=0.8,main=i*10)
}
ordispider(rt$Y,groups=site,col="grey80",alpha=0.01)
ordiellipse(rt$Y,groups= site,draw="polygon",col=colpops,label=T)

#-------------
# clustering / PCoA based on identity by state (IBS) based on single read resampling
# (for low and/or uneven coverage)

ma = as.matrix(read.table("ok.ibsMat"))
hc=hclust(as.dist(ma),"ave")
plot(hc,cex=0.5)  # this shows how similar clones are

ma=ma[goods,goods]
dimnames(ma)=list(bams[goods],bams[goods])
hc=hclust(as.dist(ma),"ave")
plot(hc,cex=0.7) # without clones

# performing PCoA and CAP
conds=data.frame(cbind(site))
pp0=capscale(ma~1)
pp=capscale(ma~site,conds)

# significance of by-site divergence
adonis(ma~site,conds)

# eigenvectors
plot(pp0$CA$eig) 

axes2plot=c(1,2)  
quartz()
library(adegenet) # for transp()
cmd=pp0
plot(cmd,choices=axes2plot,display="sites",type="n") # choices - axes to display
points(cmd,choices=axes2plot,pch=19,col=transp(colors,alpha=0.7))
#ordihull(cmd,choices= axes2plot,groups= conds$grp,draw="polygon",col=1+as.numeric(unique(as.factor(conds$grp))),label=T)
ordispider(cmd,choices= axes2plot,groups=conds$site,col="grey80")
ordiellipse(cmd,choices= axes2plot,groups= conds$site,draw="polygon",col=colpops,label=T)

# unscaled, to identify outliers
plot(cmd$CA$u[,axes2plot],pch=19,col=colors)
ordispider(cmd$CA$u[,axes2plot],groups=conds$site,col="grey80")
ordiellipse(cmd$CA$u[,axes2plot],groups= conds$site,draw="polygon",col=colpops,label=T)
identify(cmd$CA$u[,axes2plot],labels=colnames(ma),n=3,cex=0.7)




