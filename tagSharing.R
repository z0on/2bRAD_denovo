if (length(commandArgs(trailingOnly=TRUE))<1) {
	stop("
	
Tag sharing calculator / resampler
Produces distance matrices based on sharing of RAD tags between samples.

Arguments: 

infile=[filename]   Name of tab-delimited kmers table produced my mergeKmers.pl
                    rows - kmers, columns - individuals,
                    first 3 colums will be ignored (kmer id, sequence, and total count), 

nreps=3             number of replicates (random resamplings when computing similarity)

first=1             first individual to compare to all (for parallelizing)

last=[total N]      last individual to compare to all (for parallelizing)

Output: 

infile_[first]_[last]_tagshare.RData - contains two matrix objects: 
                                           strd (based on tag presence-absense)
                                           cord (based on correlations of tag counts)
        infile_[first]_[last].tagMat - tab-delimited distance matrix based on presence-absences
     infile_[first]_[last].tagCorMat - tab-delimited distance matrix based on correlations of tag counts
   
Mikhail Matz, matz@utexas.edu, 2020

")
}

infile =grep("infile=",commandArgs())
if (length(infile)==0) { stop ("specify input file (infile=filename)\n") }
infile=sub("infile=","", commandArgs()[infile])

rr=grep("nreps=",commandArgs())
fst=grep("first=",commandArgs())
lst=grep("last=",commandArgs())
if (length(rr)==0) { Nr=3 } else { Nr=as.numeric(sub("nreps=","",commandArgs()[rr])) }
if (length(fst)==0) { fst=1} else { fst=as.numeric(sub("first=","",commandArgs()[fst])) }

message(paste("reading",infile,"..."))
tgg=read.table(infile,sep="\t",header=TRUE)[,-c(1:3)] 
tgis=tgg>0
rsum=apply(tgg,2,sum)
samples=colnames(tgg)	
if (length(lst)==0) { lst=ncol(tgg)} else { lst=as.numeric(sub("last=","",commandArgs()[lst])) }

outfl=paste(infile,"_",fst,"_",lst,"_tagshare.RData",sep="")
outfl2=paste(infile,"_",fst,"_",lst,".tagMat",sep="")
outfl3=paste(infile,"_",fst,"_",lst,".tagCorMat",sep="")

# calculating tag sharing - mega long time!
message ("calculating tag share distance (prepare to wait a while)...")
stre=matrix(0,nrow=ncol(tgis),ncol=ncol(tgis))
cormat=stre
perc=seq(10,90,10);pp0=0
for (i in fst:lst) {
	if (i==ncol(tgg)) { break }
	pp=round(100*(i-fst+1)/(lst-fst),3)
#	message(i)
	for (j in (i+1):ncol(tgg)) {
#	message(paste("    ",j))
		if (sum(tgg[,i])>sum(tgg[,j])) {
			A=tgg[,i]
			B=tgg[,j]
		} else {
			A=tgg[,j]
			B=tgg[,i]
		}
		ematch=0;cors=0
		for (r in 1:Nr) {
			probs=A/sum(A)
			s1=hist(sample(c(1:length(A)),sum(B),prob=probs,replace=TRUE),breaks=c(0:length(A)),plot=F)$counts
			s2=hist(sample(c(1:length(A)),sum(B),prob=probs,replace=TRUE),breaks=c(0:length(A)),plot=F)$counts
			ematch=ematch+(1/Nr)*sum(s1*B>0)/sum(s1*s2>0)
			cors=cors+(1/Nr)*cor(s1,B)/cor(s1,s2)
		}
#	message(paste("         ematch:",ematch))
		stre[i,j]=ematch
		cormat[i,j]=cors
#	message(paste("         MATCH:",stre[i,j]))
	}
	message(paste(fst,"-",lst,":",i,"(",pp,"%)",sep="")) 
}

dimnames(stre)=dimnames(cormat)=list(samples,samples)
stre=stre+t(stre)
cormat=cormat+t(cormat)
diag(stre)=1
diag(cormat)=1
stre[stre>1]=1
cormat[cormat>1]=1
strd=1-stre
cord=1-cormat

write.table(strd,quote=FALSE, row.names=FALSE, col.names = FALSE,file=outfl2,sep="\t")
write.table(cord,quote=FALSE, row.names=FALSE, col.names = FALSE,file=outfl3,sep="\t")
save(strd,cormat,file=outfl)

