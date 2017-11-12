
if (length(commandArgs(trailingOnly=TRUE))<1) {
	stop("
	
Converts a series to ANGSD's sfs files to a dadi-format counts table
(generates counts by sampling according to probabilities for each site)
All sfs files must be generated for the same sites ( -rf option in ANGSD)

arguments: 

infiles=[filename]   list of per-population sfs files (one file per line), 
                    generated with \'realSFS mypop.saf.idx > mypop.sfs\'
prefix=[string]     prefix for _dadi.data output file

example:

sfs2dadi infiles=sfs.list prefix=mycreatures
   
Mikhail Matz, matz@utexas.edu

")
}

infile =grep("infiles=",commandArgs())
if (length(infile)==0) { stop ("specify file listing input sfs files (infiles=filename)\n") }
infile=sub("infiles=","", commandArgs()[infile])

prefix =grep("prefix=",commandArgs())
if (length(prefix)==0) { prefix="out" } else { prefix=sub("prefix=","", commandArgs()[prefix]) }

ins=as.character(read.table(infile)[,1])
fi=read.table(ins[1],sep="\t")
coord=paste(fi[,1],fi[,2],sep=".")

dadi0=matrix(nrow=length(coord),ncol=length(ins))
dimnames(dadi0)=list(coord,ins)
dadi1=dadi0
ns=2*length(ins)
for (inn in ins) {
	message (inn)
	sfs=read.table(inn,sep="\t")
	chr=sfs[,1]
	pos=sfs[,2]
	coord0=paste(chr,pos,sep=".")
	if (sum(coord0!=coord)!=0) { stop ("incompatible sfs: they must be generated for the SAME sites in all pops (use -rf option when calling angsd)") }
	sfs=sfs[,-c(1,2)]
	sfs=exp(sfs)
	cts=c()
	for(i in 1:nrow(sfs)) {
		cts=append(cts,sample(1:ncol(sfs)-1,1,prob=sfs[i,]/sum(sfs[i,])))
	}
	dadi1[,inn]=cts
	dadi0[,inn]=ncol(sfs)-cts-1
}
Allele1=rep("A",nrow(dadi1))
Allele2=rep("T",nrow(dadi1))
Gene=as.character(fi[,1])
Position=as.character(fi[,2])
NAME=rep("GAG",nrow(dadi1))
OUT=NAME
dadi=data.frame(cbind(NAME,OUT,Allele1,dadi0,Allele2,dadi1,Gene,Position))
names(dadi)=sub("\\.sfs|\\.sfs\\.1|\\.1","",names(dadi))
#head(dadi)

write.table(dadi,file=paste(prefix,"dadi.data",sep="_"),sep="\t",row.names=FALSE,quote=FALSE)

