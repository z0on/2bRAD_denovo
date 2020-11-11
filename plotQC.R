# quick script to compute percentiles and plot distributions of quality scores and depths
# cannibalized by Mikhail Matz from the original script by Matteo Fumagalli

#    Arguments:
#    prefix=[character]  : common prefix of angsd output files
#    bams=[filename]     : file listing bam files in the order they appear in the data (optional; use for angsd -pileup runs)

pref=grep("prefix=",commandArgs())
if (length(pref)==0) { stop ("specify prefix of angsd output files\n\nArguments:\nprefix=[character]  : common prefix of angsd output files\nbams=[filename]  : file listing bam files in the order they appear in the data (use for angsd -pileup runs)") }
pref=sub("prefix=","", commandArgs()[pref])

fin=pref

bams=grep("bams=",commandArgs())
if(length(bams)>0) { 
	bamlist=sub("bams=","", commandArgs()[bams])
} else {
	bamlist=system(paste("grep 'angsd.-b' ",fin,".arg | perl -pe 's/.+-b\\s(\\S+).+/$1/'",sep=""),intern=TRUE)
}

# setwd("/Users/c-monstr/Dropbox/angsd_pileup")
# pref="dd"
# bamlist="bams.combo"

bams=read.table(bamlist)[,1]

# print a list of samples in order from worst to best covered 

#cat("", file=paste(fin,".info",sep="",collapse=""))

pdf(paste(fin,".pdf",sep="",collapse=""),height=12,width=3.5)
par(mfrow=c(4,1))

cname=paste(fin,".counts", sep="")
system(paste("gunzip", paste(cname,".gz",sep="")))
cts=read.table(cname,sep="\t",header=T)
cts$X=NULL
ctss=cts>0
indcover=apply(ctss,1,sum)
indcover=indcover/ncol(ctss)
indcover=sort(indcover,decreasing=T)

## plot q-scores
qs <- read.table(paste(fin,".qs",sep="",collapse=""), head=T, stringsAsFactors=F)
qs <- cbind(qs, perc=cumsum(qs$counts/1e4)/sum(qs$counts/1e4,na.rm=T))
qs$perc=1-qs$perc
plot(perc~qscore,qs,type="s",ylab="fraction of sites remaining",main=paste(length(indcover),"sites"),xlab="base Q score cutoff")
write.table( qs, row.names=F, col.names=T, quote=F, sep="\t", file=paste(fin,".info",sep="",collapse=""), append=T)

## global depth
dep <- as.numeric(scan(paste(fin,".depthGlobal", sep="",collapse=""),what="char", quiet=T))
#barplot(height=dep, names.arg=seq(1,length(dep))-1, xlab="Global Depth", ylab="Counts")
cat("\nGlobal_depth\tpercentile\n", file=paste(fin,".info",sep="",collapse=""), append=T)
write.table( cbind(seq(1,length(dep))-1,cumsum(dep)/sum(dep)), row.names=F, col.names=F, quote=F, sep="\t", file=paste(fin,".info",sep="",collapse=""), append=T)

# sample depth
deps <- read.table(paste(fin,".depthSample", sep="",collapse=""),head=F, stringsAsFactors=F)
## per sample
depp <- matrix(NA, nrow=nrow(deps), ncol=ncol(deps))
for (i in 1:nrow(depp)) depp[i,] <- apply(X=deps[i,], FUN=sum, MAR=2)

# calculate and plot proportions of sites left after certain depth cutoff, in each sample
# this will plot only the first 'xl' bins of depth; change as appropriate
xl <- 20
depp=depp[,-1]
cumul=c()
sums=apply(depp,1,sum)
cumnum=depp[,1]
cumul=data.frame(cbind(cumul,cumnum/sums))
for (d in 2:xl){
        cumnum=(cumnum+depp[,d])
        cumul=data.frame(cbind(cumul, cumnum/sums))
}
cumul=t(cumul)
plot(1-cumul[,1],type="l",xlab="coverage cutoff", ylab="fraction of sites remaining",ylim=c(0,1),col=rgb(0,0,0,alpha=0.2))
# median proportion of sites with coverage >10, across samples
abline(v=10,lty=3)
abline(h=median(1-cumul[10,]),lty=3)
for (i in 2:ncol(cumul)) { lines(1-cumul[,i],col=rgb(0,0,0,alpha=0.2)) }

c5=1-cumul[5,]
names(c5)=bams
print("proportion of sites better than coverage of 5 for each sample:")
logq=log(c5+1e-5)
zq=(logq-mean(logq))/sd(logq)
goods=data.frame(cbind(names(zq[zq>(-3)])))
print(data.frame(cbind(sort(c5))))
write.table(data.frame(cbind(sort(c5))),file="quality.txt",quote=F,col.names=F)
write.table(goods,file="bams.qc",quote=F,row.names=F,col.names=F)

# plot number of sites left after certain genotying rate (-minInd) cutoff
sites=c(1:length(indcover))/length(indcover)
plot(sites~indcover,type="s",ylab="fraction of sites remaining",xlab="genotyping rate cutoff")
plot(sites~indcover,type="s",ylab="fraction of sites remaining",xlab="genotyping rate cutoff",log="y")

invisible(dev.off())
#setwd('~/Dropbox/Documents/ecogeno2018/mcav_classProject_2018/')

