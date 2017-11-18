if (length(commandArgs(trailingOnly=TRUE))<1) {
	stop("
	
Tag sharing calculator / resampler

arguments: 

infile=[filename]   Name of tab-delimited tag-counts table 
                    (rows - tags, columns - individuals,
                    no other columns except these, 
                    header line must be present - names of individuals)
                    generate this table in R with 
                    write.table(dataframe,sep=\"\\t\",row.names=F,quote=F,file=\"data.tab\")

boot=[1|0]          Whether to resample tags with replacement (1). Default 0.    
nreps=[integer]     Number of count-resampling replicates. At least 2 (the default)
outfile=[filename]  Name of output R data bundle containig distance matrix named strd. 
                    Default [input filename]_tagshare.RData        
   
Mikhail Matz, matz@utexas.edu

")
}

infile =grep("infile=",commandArgs())
if (length(infile)==0) { stop ("specify input file (infile=filename)\n") }
infile=sub("infile=","", commandArgs()[infile])

nreps=grep("nreps=",commandArgs())
boot=grep("boot=",commandArgs())
outfl=grep("outfile=",commandArgs())

if (length(boot)==0) { boot=0 } else { boot=sub("boot=","",commandArgs()[boot]) }
if (!(boot %in% c(0,1))) { stop ("unexpected boot setting, should be 1 or 0 (either resample or not)\n") }

if (length(nreps)==0) { nreps=2 } else { nreps=sub("nreps=","",commandArgs()[nreps]) }
if (nreps < 2) { stop ("unexpected nreps setting (number of resampling replicates), should be integer >1\n") }

if (length(outfl)==0) { outfl=paste(infile,"_tagshare.RData",sep="") } else { outfl=sub("outfile=","",commandArgs()[outfl]) }

message(paste("reading",infile,"..."))
tgg=read.table(infile,sep="\t",header=T)
if (boot==1){
	tgg=tgg[sample(1:nrow(tgg),replace=TRUE),]
}

# subsampling to same number of reads (0.75 of minimal)
sums=apply(tgg,2,sum)
minsum=min(sums) 
minsum=round(minsum*0.75,0)
message (paste("subsampling to",minsum,"tag count..."))

tgre=vector("list",nreps)
for (g in 1:ncol(tgg)) {
	 probs=tgg[,g]/sums[g]
	 for (r in 1:nreps) {
		 counts=hist(sample(c(1:nrow(tgg)),minsum,prob=probs,replace=TRUE),breaks=c(0:nrow(tgg)),plot=F)$counts
		 tgre[[r]]=data.frame(cbind(tgre[[r]],counts))
	 }
}
samples=colnames(tgg)	
tgis=vector("list",nreps)
for (r in 1:nreps){
	tgis[[r]]=tgre[[r]]>0
}

#-----------------------
# expected to-itself match for each sample

ematch=c()
nreps=length(tgis)
message ("calculating match to itself...")
ri=1;rj=2;i=2
for (i in 1:ncol(tgis[[1]])) {
	 summ=0;nn=0
	 for (ri in 1:(nreps-1)){
		 for (rj in (ri+1):nreps){
		 	nn=nn+1;
		 	summ=summ+sum(tgis[[ri]][,i] & tgis[[rj]][,i])
		 }
	 }
	 summ=summ/nn
	 ematch=append(ematch,summ)
}

# calculating tag sharing - mega long time!
message ("calculating tag share distance (prepare to wait a while)...")
nreps=length(tgis)
stre=matrix(0,nrow=ncol(tgis[[1]]),ncol=ncol(tgis[[1]]))
perc=seq(10,90,10);pp0=0
for (i in 1:(ncol(tgis[[1]])-1)) {
	pp=round(100*i/ncol(tgis[[1]]),-1)
	if (pp != pp0 & pp %in% perc) { 
		message(paste(pp,"% done")) 
		pp0=pp
	}
	for (j in (i+1):ncol(tgis[[1]])) {
		 summ=0
		 for (r in 1:nreps){
			 	summ=summ+sum(tgis[[r]][,i] & tgis[[r]][,j])
		 }
		 summ=summ/nreps	
		 stre[i,j]=summ/mean(ematch[i],ematch[j])
	}
}

dimnames(stre)=list(samples,samples)
stre=stre+t(stre)
diag(stre)=1
strd=1-stre
message(paste("min off-diag distance:",round(min(strd[lower.tri(strd)]),3),"max distance:",round(max(strd),3)))
message(paste("number of negative distances (will be converted to 0's):", sum(strd<0),"\n"))
strd[strd<0]=0
save(strd,file=outfl)

