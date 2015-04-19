# execute the chunk between ------- lines to load the function
#-----------------
plotAdmixture=function(data,npops,colors=NA,...) {
	require(RColorBrewer)
	if (is.na(colors[1])){ colors=brewer.pal(n=npops,name="Set3") }
	p=levels(data$pop)[1]
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
	midpts=barplot(t(as.matrix(tsort[1:npops])), col=colors,xlab="population", space=0,ylab="Ancestry", border=NA, xaxt="n",mgp=c(2,1,0),...)
	pops=levels(tbl$pop)
	np=0;lp=0;lpoints=c()
	abline(v=0,lwd=1)
	abline(h=1,lwd=1)
	abline(h=0,lwd=1)
	for( p in pops) {
		np0=table(tbl$pop==p)[2]
		lp0=np0/2
		lp=np+lp0
		np=np+np0
		lpoints=append(lpoints,lp)
		abline(v=np,lwd=1)
	}
	mtext(pops,side=1,line=0.5,at=lpoints)
}
#---------------------

library(RColorBrewer)

# assembling the input table
inName="admix.4.Q" # name of the input file to plot, output of ADMIXTURE run
npops=as.numeric(strsplit(inName,".",fixed=T)[[1]][length(strsplit(inName,".",fixed=T)[[1]])-1])
inds2pops="inds2pops.tab" # two-column tab-delimited file listing sample names and their population assignemnts, in the order corresponding to the ADMIXTURE output

tbl=read.table(inName)
i2p=read.table(inds2pops,header=T)
tbl=cbind(tbl,i2p)

head(tbl,20) # this is how the resulting dataset must look

cols=brewer.pal(n=npops,name="Set3")  # or use your own colors

plotAdmixture(tbl,npops=npops,colors=cols)


