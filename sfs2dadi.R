if (length(commandArgs(trailingOnly=TRUE))<1) {
	stop("
	
Converts a series to ANGSD's sfs files to a dadi-format counts table
(generates counts by sampling according to probabilities for each site)
All sfs files must be generated for the same sites ( -rf option in ANGSD)

arguments: 

infiles=[filename]   list of per-population sfs files (one file per line), 
                    generated with \'realSFS mypop.saf.idx > mypop.sfs\'
prefix=[string]     prefix for _dadi.data output file

maxsnp=1e+10	max number of SNPs to process

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

maxsnp =grep("maxsnp=",commandArgs())
if (length(maxsnp)==0) { maxsnp=1e+10 } else { maxsnp=sub("maxsnp=","", commandArgs()[maxsnp]) }
maxsnp=as.numeric(maxsnp)

message(paste("args:",infile,prefix,maxsnp))

ins=as.character(read.table(infile)[,1])
for (inn in ins) {
	message (paste ("reading",inn))
	sfs=read.table(inn,sep="\t")
	if (inn==ins[1]){
		coord=paste(sfs[,1],sfs[,2],sep=".")
		message(paste("   coord/maxsnp:",length(coord),maxsnp))
		if (length(coord)>maxsnp) { coord=coord[1:maxsnp] }
		dadi0=matrix(nrow=length(coord),ncol=length(ins))
		row.names(dadi0)=coord
		dimnames(dadi0)=list(coord,ins)
		dadi1=dadi0
		ns=2*length(ins)
	}
	sfs=sfs[1:length(coord),]
	message (paste(nrow(sfs),"SNPs",round(proc.time()[3]/60,1)))

	chr=sfs[,1]
	pos=sfs[,2]
	coord0=paste(chr,pos,sep=".")
	row.names(sfs)=coord0
	coord=coord[coord %in% coord0]
	sfs=sfs[coord,]
	dadi1=dadi1[coord,]
	dadi0=dadi0[coord,]
	if (inn!=ins[1]) { message(length(coord)," overlapping sites left")}
#	if (sum(coord0!=coord)!=0) { stop ("incompatible sfs: they must be generated for the SAME sites in all pops (use -rf option when calling angsd)") }
	sfs=sfs[,-c(1,2)]
	sfs=exp(sfs)
	cts=c()
	for(i in 1:nrow(sfs)) {
		cts=append(cts,sample(1:ncol(sfs)-1,1,prob=sfs[i,]/sum(sfs[i,])))
		if (i %in% seq(1e+4,10e+6,1e+4)) { message (paste(inn,i,round(proc.time()[3]/60,1))) }
	}
	dadi1[,inn]=cts
	dadi0[,inn]=ncol(sfs)-cts-1
}
Allele1=rep("A",nrow(dadi1))
Allele2=rep("T",nrow(dadi1))
Gene=sub("\\..+","",coord)
Position=sub(".+\\.","",coord)
#Gene=as.character(fi[,1])
#Position=as.character(fi[,2])
NAME=rep("GAG",nrow(dadi1))
OUT=NAME
dadi=data.frame(cbind(NAME,OUT,Allele1,dadi0,Allele2,dadi1,Gene,Position))
names(dadi)=sub("\\.sfs|\\.sfs\\.1|\\.1","",names(dadi))

write.table(dadi,file=paste(prefix,"dadi.data",sep="_"),sep="\t",row.names=FALSE,quote=FALSE)

