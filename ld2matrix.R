fin = commandArgs(T)

message("reading input...")
ald=read.table(fin,sep="\t") #[,c(1,2,7)]
message("...done")
names(ald)=c("s1","s2","rEM")
sites=sort(unique(c(as.character(ald$s1),as.character(ald$s2))))
nsites=length(sites)
ldmat=matrix(-1,nrow=nsites,ncol=nsites)
dimnames(ldmat)=list(sites,sites)
message("analyzing sites...")
for (i in 1:(nsites-1)) {
	message(paste(i,":",sites[i]))
	sub=subset(ald,s1==sites[i] | s2==sites[i])
	swi=which(sub$s2==sites[i])
	sub$s2[swi]=sub$s1[swi]
#	print(sub[1:5,])
	for (j in which(ldmat[,i]==-1)) {
		if (i==j) { next }
#		message("    ",j,sub$s2[j])
	 	ldmat[sites[i],sub$s2[j]]=ldmat[sub$s2[j],sites[i]]=sub[j,"rEM"]
	}
}
diag(ldmat)=1
save(ldmat,file=paste(fin,"_matrix.RData",sep="",collapse=""))