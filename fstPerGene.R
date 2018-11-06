# reading gene regions
genes=read.table("genes.txt")
names(genes)=c("contig","start","end","gene")

# adding gene annotations, in the same order as genes.txt (skip if you don't have such file):
gnames=read.table("gnames.txt",sep="\t")
names(gnames)=c("gene","protein")
genes=merge(genes,gnames,by="gene",all.x=T)
genes$protein=as.character(genes$protein)
genes$protein[is.na(genes$protein)]="unknown"
nrow(genes[genes$protein!="unknown",])

fst01=read.table("p01.fst")
names(fst01)=c("contig","pos","a","b")

fst01$contig=as.character(fst01$contig)
fst34$contig=as.character(fst34$contig)

# removing zero-only (invariant) bases
fst01[,3:4]=round(fst01[,3:4],3)
ch01=apply(fst01[,3:4],1,sum)
chh=(ch01>0)
table(chh)
 # FALSE   TRUE 
# 932541  41886 

fst01=fst01[chh,]
head(fst01)

# Fst density
jitter=0
#jitter=rnorm(nrow(fst01),0,0.0)
f1=fst01$a/(fst01$b+1e-3)+jitter
plot(density(f1),ylim=c(0,7),xlab="Fst",main="",mgp=c(2.3,1,0),bty="n")

# computing weighted Fst per gene
i=1;gfst01=c();ns01=0
pb=txtProgressBar(0,nrow(genes))
for (i in 1:nrow(genes)) {
	setTxtProgressBar(pb,i)
	sub=subset(fst01,contig==genes$contig[i] & pos>=genes$start[i] & pos<=genes$end[i])
	if (is.null(sub[1,1]) | sum(sub$b)==0) { 
		gfst01=append(gfst01,NA)
	} else {
		gfst01=append(gfst01,sum(sub$a)/sum(sub$b))
		ns01=ns01+nrow(sub)
	}
}

# density plots of per-gene Fst
plot(density(na.omit(gfst01)))

# adding results to genes table, saving
genes$fst01=gfst01
head(genes[genes$protein!="unknown",])
save(genes,file="PerGeneFst.RData")

# saving table of per-gene Fst for GO_MWU analysis
f01=genes[,c("gene","fst01")]
f01=na.omit(f01)
nrow(f01)
write.csv(f01,file="f01.csv",quote=F,row.names=F)

