infile=commandArgs(T)
if (length(infile)<2) {
	message("thinner_dadi: resamples sites from dadi-formatted file (first arg)
	until minimum specified distance (arg2) is approached.
	usage example: Rscript thinner_dadi.R dadi.data 5000\n")
	stop("specify input file and distance\n")
}

#infile=c("OK_dadi.data","5000")

sfs=read.table(infile[1],header=T)
di=as.numeric(infile[2])

resfs=c();g="chr1"
for (g in levels(sfs$Gene)) {
	message(g)
	ss=subset(sfs,Gene==g)
	poss=c(as.numeric(sample(ss$Position,1)))
	lesser=0;
	while (lesser <500) {
		pos=as.numeric(sample(ss$Position,1));p=1;close=0
		for (p in 1:length(poss)) {
			if (abs(poss[p]-pos)<di) {
				close=1
			}
		}
		if (close==1) { 
				lesser=lesser+1
		} else {
				poss=append(poss,pos)
				lesser=0
		}
	}
	ss=subset(ss,Position %in% poss)
	resfs=data.frame(rbind(resfs,ss))
}
message(paste(nrow(sfs),"sites total"))
message(paste(nrow(resfs),"selected"))
write.table(resfs,file=paste("thin",di,infile[1],sep="_"),sep="\t",quote=F,row.names=F)

