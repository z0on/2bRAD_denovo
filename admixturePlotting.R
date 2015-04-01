library(RColorBrewer)
#----------------------------------------
popline=function(poplabels,midpoints) {	
	p0=poplabels[1];pi=0;popIDpts=c();tline=-0.4
	for (p in 1:length(poplabels)) {
		 if (poplabels[p] != p0 ) { 
		 		if (tline==-0.4){tline=-0.1}else{tline=-0.4}
			 	popIDpts=append(popIDpts,p-ceiling(pi/2)-1)
			 	p0=poplabels[p]
			 	pi=0
		 } else { pi=pi+1 } 
		 mtext("-",side=1,line=tline,at=midpoints[p])
	}
	popIDpts=append(popIDpts,p-floor(pi/2))
	for (p in 1:length(unique(poplabels))){
		  mtext(unique(poplabels)[p],side=1,line=0.5,at=midpoints[popIDpts[p]])
	}
}#----------------------------------------
# plink --vcf denov.vcf --make-bed --out denov  # (must have plink 1.09)
# grep "#CHROM" denov.vcf | perl -pe 's/\t/\n/g' | grep [0-9] > inds.list
# admixture denov.bed 5

tbl=read.table("denov.5.Q")
inds=read.table("inds.list")
tbl$ind=inds$V1
tbl$pop=gsub("[0-9]","",as.character(inds$V1))
tbl$pop=factor(tbl$pop,levels=c("Ni","S","O","M","K"))
npops=5
p=levels(tbl$pop)[1]
tsort=c()
for(p in levels(tbl$pop)){
	s=subset(tbl,pop==p)
	probs=apply(s[,1:npops],2,sum)
	ord=order(probs,decreasing=T)
	if (npops>=5){so=s[order(s[,ord[1]],s[,ord[2]],s[,ord[3]],s[,ord[4]],s[,ord[5]],decreasing=T),]} 
	if (npops==4){ so=s[order(s[,ord[1]],s[,ord[2]],s[,ord[3]],s[,ord[4]],decreasing=T),]} 
	if (npops==3){ so=s[order(s[,ord[1]],s[,ord[2]],s[,ord[3]],decreasing=T),]} 
	if (npops==2) { so=s[order(s[,ord[1]],s[,ord[2]],decreasing=T),] }
	tsort=data.frame(rbind(tsort,so))
}

cols=brewer.pal(n=5,name="Accent")
midpts=barplot(t(as.matrix(tsort[1:npops])), col=cols,xlab="population", space=0,ylab="Ancestry", border=NA, xaxt="n",mgp=c(2.3,1,0))
popline(tsort$pop,midpts)
