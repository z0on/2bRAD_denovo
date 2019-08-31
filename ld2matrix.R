fin = commandArgs(T)
fin="LD_1.LD.slim"
message("reading input...")
ald=read.table(fin,sep="\t",stringsAsFactors=FALSE) #[,c(1,2,7)]
message("...done")
names(ald)=c("s1","s2","rEM")
sites=sort(unique(c(ald$s1,ald$s2)))
nsites=length(sites)
ldmat=matrix(-1,nrow=nsites,ncol=nsites)
dimnames(ldmat)=list(sites,sites)
message("analyzing sites...")
for (i in 1:nsites) {
	message(paste(i,":",sites[i]))
	sub=subset(ald,s1==sites[i] | s2==sites[i])
	nrow(sub)
	swi=which(sub$s2==sites[i])
	sub$s2[swi]=sub$s1[swi]
	row.names(sub)=sub$s2
	for (j in which(ldmat[,i]==-1)) {
		if (i==j) { next }
	 	ldmat[sites[i],sites[j]]=ldmat[sites[j],sites[i]]=sub[sites[j],"rEM"]
	}
}
diag(ldmat)=1
save(ldmat,file=paste(fin,"_matrix.RData",sep="",collapse=""))