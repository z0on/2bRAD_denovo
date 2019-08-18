if (length(commandArgs(trailingOnly=TRUE))<1) {
	stop("
	
Tag sharing calculator, version 2

arguments: 

infile=[filename]   Name of RData file containing tag-counts table 
                    (rows - tags, columns - individuals,
                    no other columns except these, 
                    header line must be present - names of individuals)

first=[1]             first individual to compare to all (for parallelizing)

last=[total N]      last individual to compare to all (for parallelizing)

nreps=[1]           number of resampling replicates 
                    when computing expected distance under null

outfile=[filename]  Name of output R data bundle containig tag-sharing matrix named stre. 
                    Default [input filename]_[first]_[last]_tagshare.RData
                    (convert this matrix to distances: 1-stre, diag=0)
                    
Mikhail Matz, 2019
matz@utexas.edu

")
}

infile =grep("infile=",commandArgs())
if (length(infile)==0) { stop ("specify input file (infile=filename)\n") }
infile=sub("infile=","", commandArgs()[infile])

outfl=grep("outfile=",commandArgs())
rr=grep("nreps=",commandArgs())
fst=grep("first=",commandArgs())
lst=grep("last=",commandArgs())
if (length(rr)==0) { Nr=1 } else { Nr=as.numeric(sub("nreps=","",commandArgs()[rr])) }
if (length(fst)==0) { fst=1} else { fst=as.numeric(sub("first=","",commandArgs()[fst])) }

message(paste("reading",infile,"..."))
load(infile) 
tgis=tgg>0
rsum=apply(tgg,2,sum)
samples=colnames(tgg)	
if (length(lst)==0) { lst=ncol(tgg)} else { lst=as.numeric(sub("last=","",commandArgs()[lst])) }
if (length(outfl)==0) { outfl=paste(infile,"_",fst,"_",lst,"_tagshare.RData",sep="") } else { outfl=sub("outfile=","",commandArgs()[outfl]) }

# calculating tag sharing - mega long time!
message ("calculating tag share distance (prepare to wait a while)...")
stre=matrix(0,nrow=ncol(tgis),ncol=ncol(tgis))
perc=seq(10,90,10);pp0=0
for (i in fst:lst) {
	if (i==ncol(tgg)) { break }
	pp=round(100*(i-fst+1)/(lst-fst),3)
#	message(i)
	for (j in (i+1):ncol(tgg)) {
#	message(paste("    ",j))
		A=tgg[,i]+tgg[,j]
		ematch=0
		for (r in 1:Nr) {
			probs=A/sum(A)
#	message(paste("      rep",r,"i"))
			si=hist(sample(c(1:length(A)),sum(tgg[,i]),prob=probs,replace=TRUE),breaks=c(0:length(A)),plot=F)$counts
#	message(paste("      rep",r,"j"))
			sj=hist(sample(c(1:length(A)),sum(tgg[,j]),prob=probs,replace=TRUE),breaks=c(0:length(A)),plot=F)$counts
#	message(paste("      rep",r,"done"))
			ematch=ematch+sum(si>0 & sj>0)/Nr
		}
#	message(paste("         ematch:",ematch))
		stre[i,j]=sum(tgis[,i] & tgis[,j])/ematch
#	message(paste("         MATCH:",stre[i,j]))
	}
	message(paste(fst,"-",lst,":",i,"(",pp,"%)",sep="")) 
}

dimnames(stre)=list(samples,samples)
stre=stre+t(stre)
# diag(stre)=1
# strd=1-stre
# message(paste("min off-diag distance:",round(min(strd[lower.tri(strd)]),3),"max distance:",round(max(strd),3)))
# message(paste("number of negative distances (will be converted to 0's):", sum(strd<0),"\n"))
# strd[strd<0]=0
save(stre,file=outfl)

