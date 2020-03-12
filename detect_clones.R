# detects clones by looking at hierarchical clustering of samples based on IBS from ANGSD
# clones are samples that cluster at very low height, < cutoff (typically 0.15-0.18)
# run this once to make hctree.pdf, check where you want to put cutoff

args=commandArgs(T)
bamfile=args[1]
ibs=args[2]
cutoff=args[3]

bams=read.table(bamfile)[,1] # list of bam files
ma = as.matrix(read.table(ibs))
dimnames(ma)=list(bams,bams)
hc=hclust(as.dist(ma),"ave")
pdf("hctree.pdf",height=8, width=15)
plot(hc,cex=0.7) 
abline(h=cutoff, col="red",lty=3)
dev.off()
cuts=cutree(hc,h=cutoff)
nrbams=c()
for (i in unique(cuts)) { 
	cc=cuts[cuts==i]
	nrbams=append(nrbams,sample(names(cc),1)) 
}
nrn=length(nrbams)
nn=length(bams)
if (nrn<nn) { 
	message("detect_clones: highly similar samples detected at cutoff ",cutoff)		
	message("check hctree.pdf")	
}
message("retained ", nrn," out of ",nn," bams")	
write(nrbams,file="bams.nr")
