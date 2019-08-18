# reading tags table
tg=read.table("allTags.uniq",sep="\t",header=T)

# creating a tag presence-absence matrix
tgg=tg[,-c(1:4)]
tgis=tgg>0

# individuals per tag:
indcover=apply(tgis,1,sum)
hist(log(indcover,10))

# tags per individual:
tperi=apply(tgis,2,sum)
hist(log(tperi,10),breaks=50)

# reads per individual:
rperi=apply(tgg,2,sum)
hist(log(rperi,10),breaks=50)

plot(tperi~rperi)
plot(tperi~rperi,log="xy")

save(tgg,file="tgAll.RData")

#---------------- calculating chunks for parlellizing calculations

x=1:ncol(tgg)
nc=ncol(tgg)-x
ncomp=sum(nc)
cum=c();thrds=36;chunks=list();start=1
th=1
for (i in x) {
	cum[i]=sum(nc[1:i])
	if(cum[i]>(ncomp/thrds)*th) {
		chunks[[th]]=data.frame(start=start,end=i)
		th=th+1
		start=i+1
	}
}
chunks[[th]]=data.frame(start=start,end=ncol(tgg))

# printing commands to run:
for (i in 1:th) {
message(paste("Rscript TagSharing_auto.R infile=tgAll.RData nreps=3 first=",chunks[[i]]$start," last=",chunks[[i]]$end,sep=""))
}
# copy the printed commands, scp TagSharing_v2.R and tgAll.tab to a cluster, run them in parallel (note: if the tag table is large, each task might require 2-4G of RAM)
# scp all the resulting *tagshare.RData files back to laptop

#--------------- assembling results

ts=matrix(0,nrow=ncol(tgis),ncol=ncol(tgis))
samples=colnames(tgis)
dimnames(ts)=list(samples,samples)
i=2
library(pheatmap)
for (i in 1:th) {
	fname=paste("tgALL.RData_",chunks[[i]]$start,"_",chunks[[i]]$end,"_tagshare.RData",sep="")
	message(fname)
	load(fname)
	ts=ts+stre
#	ts[1:20,1:20]
#	pheatmap(ts,cluster_rows=F,cluster_cols=F)
}
diag(ts)=1
strd=1-ts
message(paste("min off-diag distance:",round(min(strd[lower.tri(strd)]),3),"max distance:",round(max(strd),3)))
message(paste("number of negative distances (will be converted to 0's):", sum(strd<0),"\n"))
strd[strd<0]=0

#---- hierarchical clustering, to display where clones are

hc=hclust(as.dist(strd),"ave")
plot(hc,cex=0.6)

#---- coloring tree by read cover

require(sparcl)

cover=rep("red",nrow(strd))
cover[rperi>quantile(rperi,0.1)]="gold"
cover[rperi>quantile(rperi,0.2)]="black"
#pdf("colorclones_byCover.pdf",width=80,height=10)
ColorDendrogram(hc, y = cover,branchlength=0.2,labels=colnames(strd),cex=0.6)
#dev.off()

#------------
# pcoa and "constrained analysis of proximities" (cap)

# loading individual to population correspondences
i2p=read.table("inds2pops",sep="\t")
row.names(i2p)=i2p[,1]
site=i2p[,2]

# setting color scheme
palette(rainbow(length(unique(site))))
colors=as.numeric(as.factor(site))
colpops=as.numeric(as.factor(sort(unique(site))))

# calculating ordination
library(vegan) 
conds=data.frame(cbind(site))
# unconcstrained (=PCoA)
pp0=capscale(as.dist(strd)~1)
# constrained (displaying variation of interest first)
pp=capscale(as.dist(strd)~site,conds)

# testing if site is significant
adonis(as.dist(strd)~site,conds)

# how many intereting PCs we have in unconstrained PCoA
plot(pp0$CA$eig)

axes2plot=c(1,2)  # try 3,2 too
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
identify(cmd$CA$u[,axes2plot],labels=colnames(stnr2),n=3,cex=0.7)

#-------------
# t-SNE:  machine learning to identify groups of samples

library(Rtsne)
library(vegan)
library(adegenet)
quartz()

# perplexity: expected number of neighbors. Set to 0.5x of N samples/pop
perp=10
# execute to the end, watch the learning process on the plot
rt = Rtsne(as.dist(strd), perplexity=perp,max_iter=2,is_distance=T)
for (i in 1:250){
	rt = Rtsne(as.dist(strd), perplexity=perp,max_iter=10,Y_init=rt$Y,is_distance=T)
	plot(rt$Y,col=colors, pch=16,cex=0.8,main=i*10)
}
ordispider(rt$Y,groups=site,col="grey80",alpha=0.01)
ordiellipse(rt$Y,groups= site,draw="polygon",col=colpops,label=T)



