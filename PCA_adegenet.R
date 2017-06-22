#-----------------
# reading VCF into genlight format
install.packages("vcfR")
install.packages("adegenet")
install.packages("vegan")

# replace with path where your vcf file is
setwd('~/Documents/mega2017/2bRAD/')
library(vcfR)
library(adegenet)
library(vegan) # for ordispider, ordihull, ordiellipse

#----------------------
# reading vcf file 

# replace vcf name in next 2 lines with yours
gl=vcfR2genlight(read.vcfR("thinMaxaf.vcf"))

# (assumes that first letter in sample name identifies population, if not, make your own inds2pops table)
system("grep '#CHROM' thinMaxaf.vcf | perl -pe 's/\\t/\\n/g' | tail -n +10 | perl -pe 's/^(\\S)(.+)/$1$2\t$1/' > inds2pops")

# or:
pops=read.table("inds2pops",sep="\t")
pop(gl)=pops$V2

#-------------
# finding samples with too much missing data

nas=c()
inds=c()
for (i in 1:nrow(gl)){
	snps=gl[i]$gen
	nas=append(nas,length(snps[[1]]$NA.posi))
	inds=append(inds,gl$ind.names[i])
}
nas=nas/gl$n.loc
hist(nas,breaks=20)

# conisder removing outliers:
# worst genotyped:
gl@ind.names[which(nas==max(nas))]

# worse than a set cutoff for missing data:
cutoff=0.1
gl@ind.names[which(nas>cutoff)]
# K3, O3

# restricting plot to those passing missing data cutoff
gl=gl[nas<cutoff] 

#-----------------------
# plotting 500 randomly chosen SNPs to check for weird patterns

gl.s=gl[,sample(1:ncol(gl),500)]
plot(gl.s,col=c("grey70","coral","purple"),legend=F)

#-----------------------
# making PCA

pca=glPca(gl,nf=3,parallel=F)

plot(pca$scores[,1:2],col=transp(as.numeric(as.factor(pop(gl))),0.3),pch=19)
ordispider(pca$scores[,1:2],pop(gl),col=as.numeric(as.factor(pop(gl))),label=T)
#ordiellipse(pca$scores[,1:2],pops$V2,label=T,draw="polygon",col="grey90")
#ordihull(pca$scores[,1:2],pops$V2,label=T,draw="polygon",col="grey90",cex=2)

# manually identifying outliers: click on outlier points in the plot.
# adjust n=3 parameter to identify more than 3 points
outliers=identify(pca$scores[,1:2],labels=gl@ind.names,n=3,cex=0.8)

# re-making PCA without outliers
gl2=gl[-outliers]
pca2=glPca(gl2,nf=3,parallel=F)
colors=as.numeric(as.numeric(as.factor(levels(pop(gl))))) # or use your own colors for populations
s.class(pca2$scores[],pop(gl2),col=transp(colors,0.5),cstar=1,cellipse=1,clabel=1,axesell=F,grid=F,cpoint=2)


