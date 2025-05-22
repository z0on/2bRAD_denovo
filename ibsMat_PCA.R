library(ggplot2)
library(vegan)

qu=read.table("quality.txt") # output of plotQC.R 
bams=qu$V1
qual=qu$V2

ma = as.matrix(read.table("mydata.ibsMat")) # output of angsd -makeIBS 1
dimnames(ma)=list(bams,bams)

# hierarchical clustering to detect wrong collectons, clones, and closely related individuals
 
hc=hclust(as.dist(ma),"ave")
plot(hc,cex=0.5) 

# assuming clones and close relatives are removed:

# performing PCoA
pp0=capscale(ma~1)

# plot eigenvalues (% of variation explained)
plot(100*pp0$CA$eig/sum(pp0$CA$eig)) 

# plotting PCs
axes2plot=c(1,2)
scores=data.frame(scores(pp0,scaling="sites",display="sites",choices=axes2plot))

# colored by depth-quaity
ggplot(scores,aes(MDS1,MDS2,color=qual))+geom_point()+theme_minimal()



