inf=commandArgs(T)

require(vegan)
ma = as.matrix(read.table(inf))

bams=scan("bams.nr",what="character") # list of bam files
bams=sub(".bam","",bams)
qc=read.table("qualRanks",skip=2,header=F)
names(qc)=c("srr","q")
qc$srr=sub(".bam","",qc$srr)
row.names(qc)=qc$srr
qc=qc[qc$srr %in% bams,]
qc=qc[bams,]

pp0=capscale(ma~qc$q)

# selecting [non-outlier] samples within 0.05-0.95 quantile range of PC1
pc1=pp0$CA$u[,1]
goods=which(pc1<quantile(pc1,0.95) & pc1>quantile(pc1,0.05))
qc=qc[goods,]

pp1=capscale(ma[goods,goods]~qc$q)
pdf(file=paste(inf,"_pcoa_eigens.pdf",sep=""),height=4.8,width=8)
par(mfrow=c(1,2))
plot(MDS2~MDS1,pp1$CA$u,asp=1)
eigsn=pp1$CA$eig/sum(pp1$CA$eig)
plot(eigsn)
lines(bstick(nrow(ma)-1),col="red")
dev.off()
print(round(eigsn[1]/bstick(nrow(ma)-1)[1],2))

