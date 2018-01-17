
# quick script to compute percentiles and plot distributions of quality scores and depths

fin <- commandArgs(T)

cat("", file=paste(fin,".info",sep="",collapse=""))

pdf(paste(fin,".pdf",sep="",collapse=""))
par(mfrow=c(2,1))

## barplot q-scores
qs <- read.table(paste(fin,".qs",sep="",collapse=""), head=T, stringsAsFactors=F)
barplot(height=qs$counts, names.arg=qs$qscore, xlab="Q-score", ylab="Counts")
qs <- cbind(qs, perc=cumsum(qs$counts/1e4)/sum(qs$counts/1e4,na.rm=T))
write.table( qs, row.names=F, col.names=T, quote=F, sep="\t", file=paste(fin,".info",sep="",collapse=""), append=T)

## global depth
dep <- as.numeric(scan(paste(fin,".depthGlobal", sep="",collapse=""),what="char", quiet=T))
#barplot(height=dep, names.arg=seq(1,length(dep))-1, xlab="Global Depth", ylab="Counts")
cat("\nGlobal_depth\tpercentile\n", file=paste(fin,".info",sep="",collapse=""), append=T)
write.table( cbind(seq(1,length(dep))-1,cumsum(dep)/sum(dep)), row.names=F, col.names=F, quote=F, sep="\t", file=paste(fin,".info",sep="",collapse=""), append=T)

## sample depth
deps <- read.table(paste(fin,".depthSample", sep="",collapse=""),head=F, stringsAsFactors=F)
## per sample
depp <- matrix(NA, nrow=nrow(deps), ncol=ncol(deps))
for (i in 1:nrow(depp)) depp[i,] <- apply(X=deps[i,], FUN=sum, MAR=2)

# this will plot only the first 'xl' bins of depth; change as appropriate
xl <- 20
depp=depp[,-1]
#barplot(depp[,1:xl], names.arg=c(seq(1,ncol(deps))-1)[1:xl], xlab="Sample Depth", ylab="Counts", beside=T)

cumul=c()
sums=apply(depp,1,sum)
cumnum=depp[,1]
cumul=data.frame(cbind(cumul,cumnum/sums))
for (d in 2:xl){
        cumnum=(cumnum+depp[,d])
        cumul=data.frame(cbind(cumul, cumnum/sums))
}
cumul=t(cumul)
plot(1-cumul[,1],type="l",xlab="sample depth", ylab="fraction of sites",ylim=c(0,1),col=rgb(0,0,0,alpha=0.2))
abline(v=10,lty=3)
abline(h=median(1-cumul[10,]),lty=3)
for (i in 2:ncol(cumul)) { lines(1-cumul[,i],col=rgb(0,0,0,alpha=0.2)) }
invisible(dev.off())


