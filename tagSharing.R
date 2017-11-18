# reading tags table
tg=read.table("5to50.uniq",sep="\t",header=T)
# histogram of counts
hist(log(tg$count,10),breaks=150)

# creating a tag presence-absence matrix
tgis=tg
row.names(tgis)=tgis[,1]
tgis=tgis[,-c(1:4)]
tgis=tgis>0
head(tgis)

# individuals per tag:
indcover=apply(tgis,1,sum)

# tags per individual:
tperi=apply(tgis,2,sum)
hist(tperi,breaks=20)
# recording outliers (too few tags) - change 40000 in the next line to whatever number makes sense for you as the low cutoff 
outliers=names(tperi)[tperi<40000]

#------------
# removing outliers and filtering

tg1=tg[,!(colnames(tg) %in% outliers)]
tgis1=tgis[,!(colnames(tgis) %in% outliers)]

indcover=apply(tgis1,1,sum)
# retaining only tags shared among at least 5 and at most 20 individuals - change 20 to 1.5x number of samples per population in your study
goods=which(indcover>=5 & indcover<=20 & tg1$count>10 & tg1$count<3000)
length(goods) 

tgis2=tgis1[goods,]
tg2=tg1[goods,]

save(tgis2,tg2,file="cleanedRaw_5to20.RData")
#----------------
# calculating tag sharing
load("cleanedRaw_5to20.RData")
write.table(tg2[,5:ncol(tg2)],sep="\t",row.names=F,quote=F,file="tgAll_5to20.tab")

system('Rscript TagSharing_auto.R infile=tgAll_5to20.tab boot=0 nreps=3')
load("tgAll_5to20.tab_tagshare.RData")

# hierarchical clustering, to display where clones are
hc=hclust(as.dist(strd),"ave")
plot(hc,cex=0.7)

#-------------
# re-scaling based on clones (skip this if you don't have genotyping replicates)

load("tgAll_5to20.tab_tagshare.RData")

# reading table of pairs of replicates 
clonepairs=read.table("clonepairs.tab",sep="\t")
repsa= clonepairs[,1]
repsb= clonepairs[,2]

# computing median similarity between clones
cd=c()
stre=1-strd
for (i in 1:length(reps)) {
	cd=append(cd,stre[repsa[i],repsb[i]])
}
cloneSim=median(cd)

# rescaling to clones similarity, equating all clone pairs to similarity of 1
stre2=stre/cloneSim
diag(stre2)=1
for (i in 1:length(reps)) {
	stre2[repsa[i],repsb[i]]=1
	stre2[repsb[i],repsa[i]]=1
}
strd2=1-stre2

# rescaled clustering tree:
hc=hclust(as.dist(strd2),"ave")
plot(hc,cex=0.5)

#-------------
# removing extra clones 

Goods=which(!(row.names(strd2) %in% repsb))
stnr2=strd2[Goods,Goods]

# how does clustering tree look without clonemates?
hc=hclust(as.dist(stnr2),"ave")
plot(hc,cex=0.5)

save(stnr2,file="5to20_rescaled_NR.RData")

#------------
# pcoa and "constrained analysis of proximities" (cap)

ll=load("5to20_rescaled_NR.RData")

# loading individual to population correspondences
i2p=read.table("inds2pops",sep="\t")
row.names(i2p)=i2p[,1]
i2p=i2p[Goods,]
site=i2p[,2]

# setting color scheme
palette(rainbow(length(unique(site))))
colors=as.numeric(as.factor(site))
colpops=as.numeric(as.factor(sort(unique(site))))

# calculating ordination
library(vegan) 
conds=data.frame(cbind(site))
# unconcstrained (=PCoA)
pp0=capscale(as.dist(stnr2)~1)
# constrained (displaying variation of interest first)
pp=capscale(as.dist(stnr2)~site,conds)

# testing if site is significant
adonis(as.dist(stnr2)~site,conds)
# site       2    0.7688 0.38438  3.7254 0.1156  0.001 ***
# Residuals 57    5.8812 0.10318         0.8844           

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
rt = Rtsne(as.dist(stnr2), perplexity=perp,max_iter=2,is_distance=T)
for (i in 1:250){
	rt = Rtsne(as.dist(stnr2), perplexity=perp,max_iter=10,Y_init=rt$Y,is_distance=T)
	plot(rt$Y,col=colors, pch=16,cex=0.8,main=i*10)
}
ordispider(rt$Y,groups=site,col="grey80",alpha=0.01)
ordiellipse(rt$Y,groups= site,draw="polygon",col=colpops,label=T)



